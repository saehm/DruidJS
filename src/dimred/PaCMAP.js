import { wasmKnnBlock } from "../matrix/distance_matrix.wasm.js";
import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";
import { PCA } from "./PCA.js";

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {KNN} from "../knn/KNN.js" */
/** @import {ParametersPaCMAP} from "./index.js" */

/** Rows per block of the exact neighbour search, chosen so one block stays cache-resident. */
const KNN_BLOCK_ROWS = 512;

/** Shared empty rejection list, so sampling with nothing to reject allocates no subarray. */
const EMPTY = new Int32Array(0);

/**
 * Pairwise Controlled Manifold Approximation Projection (PaCMAP)
 *
 * A dimensionality reduction technique that uses three types of point pairs — nearest neighbor
 * (NN), mid-near (MN), and further (FP) pairs — with a dynamic three-phase weight schedule and
 * Adam optimization to preserve both local and global structure. Neighbours are scaled by a local
 * density estimate before selection, which is what separates the NN set from a plain KNN graph.
 *
 * @class
 * @template {InputType} T
 * @template {ParametersPaCMAP} [P=ParametersPaCMAP] - Widened by {@link LocalMAP}, which adds a parameter.
 * @extends DR<T, P>
 * @category Dimensionality Reduction
 * @see {@link https://arxiv.org/abs/2012.04456|PaCMAP Paper}
 * @see {@link https://github.com/YingfanWang/PaCMAP|PaCMAP GitHub}
 * @see {@link UMAP} for a related graph-based technique
 * @see {@link LocalMAP} for the local-refinement variant
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const X = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]];
 * const pacmap = new druid.PaCMAP(X, {
 *     n_neighbors: 10,
 *     MN_ratio: 0.5,
 *     FP_ratio: 2.0,
 *     seed: 42
 * });
 *
 * const Y = pacmap.transform(); // 450 iterations (default)
 * // [[x1, y1], [x2, y2], [x3, y3]]
 */
export class PaCMAP extends DR {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<P>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters) {
        super(
            X,
            /** @type {P} */ ({
                n_neighbors: 10,
                MN_ratio: 0.5,
                FP_ratio: 2.0,
                d: 2,
                metric: euclidean,
                lr: 1.0,
                num_iters: [100, 100, 250],
                knn: null,
                apply_pca: true,
                seed: 1212,
            }),
            parameters,
        );
        [this._N, this._D] = this.X.shape;
        const n_neighbors = /** @type {number} */ (this.parameter("n_neighbors"));
        if (n_neighbors >= this._N) {
            throw new Error(
                `Parameter n_neighbors (=${n_neighbors}) needs to be smaller than dataset size (N=${this._N})!`,
            );
        }
        this._iter = 0;
    }

    /**
     * Finds the `n_candidates` nearest neighbors of every point.
     *
     * The reference over-fetches 50 beyond `n_neighbors` so that the density rescaling in
     * {@link _scaled_neighbors} has something to re-rank, which makes this the most expensive part
     * of a run — asking an approximate index for 60 neighbors instead of 10 costs it far more than
     * the extra 50 suggests. The default path is therefore an exact blocked search whose selection
     * happens inside the WASM kernel; it beats a tree at these widths and matches the exact search
     * the reference performs. A `knn` index is used instead when one is supplied.
     *
     * @protected
     * @param {Matrix} X - The matrix to search, already reduced if PCA was applied.
     * @param {number} n_candidates
     * @returns {{ dist: Float64Array; idx: Int32Array }} Row-major `N` ⨯ `n_candidates`, ascending.
     */
    _find_candidates(X, n_candidates) {
        const N = this._N;
        const D = X.shape[1];
        const knn = /** @type {KNN<number[] | Float64Array, any> | null} */ (this.parameter("knn"));
        const dist = new Float64Array(N * n_candidates);
        const idx = new Int32Array(N * n_candidates);

        if (knn) {
            for (let i = 0; i < N; ++i) {
                const found = knn.search_by_index(i, n_candidates);
                if (found.length < n_candidates) {
                    throw new Error(
                        `The KNN index returned only ${found.length} of the ${n_candidates} requested neighbors for element ${i}. ` +
                            `Lower "n_neighbors", or use an index that recalls more candidates.`,
                    );
                }
                for (let c = 0; c < n_candidates; ++c) {
                    dist[i * n_candidates + c] = found[c].distance;
                    idx[i * n_candidates + c] = found[c].index;
                }
            }
            return { dist, idx };
        }

        // The kernel is euclidean-only, so anything else falls back to the JS search below.
        const metric = /** @type {Metric} */ (this.parameter("metric"));
        if (metric === euclidean) {
            const row_len = 2 * n_candidates;
            const out = new Float64Array(N * row_len);
            let ok = true;
            for (let start = 0; start < N && ok; start += KNN_BLOCK_ROWS) {
                const end = Math.min(start + KNN_BLOCK_ROWS, N);
                ok = wasmKnnBlock(X.values, N, D, n_candidates, out, start, end);
            }
            if (ok) {
                for (let i = 0; i < N; ++i) {
                    for (let c = 0; c < n_candidates; ++c) {
                        dist[i * n_candidates + c] = out[i * row_len + c];
                        idx[i * n_candidates + c] = out[i * row_len + n_candidates + c];
                    }
                }
                return { dist, idx };
            }
        }

        return this._find_candidates_js(X, n_candidates);
    }

    /**
     * Exact neighbour search in JS, for a non-euclidean metric or when WASM is unavailable.
     *
     * @protected
     * @param {Matrix} X
     * @param {number} n_candidates
     * @returns {{ dist: Float64Array; idx: Int32Array }}
     */
    _find_candidates_js(X, n_candidates) {
        const N = this._N;
        const metric = /** @type {Metric} */ (this.parameter("metric"));
        const dist = new Float64Array(N * n_candidates);
        const idx = new Int32Array(N * n_candidates);
        // Bounded max-heap, so the scan stays O(N log k) per row instead of sorting all of N.
        const hd = new Float64Array(n_candidates);
        const hi = new Int32Array(n_candidates);

        for (let i = 0; i < N; ++i) {
            const x_i = X.row(i);
            let size = 0;
            for (let j = 0; j < N; ++j) {
                if (j === i) continue;
                const v = metric(x_i, X.row(j));
                if (size < n_candidates) {
                    hd[size] = v;
                    hi[size] = j;
                    let c = size++;
                    while (c > 0) {
                        const p = (c - 1) >> 1;
                        if (hd[p] >= hd[c]) break;
                        const t = hd[p];
                        hd[p] = hd[c];
                        hd[c] = t;
                        const u = hi[p];
                        hi[p] = hi[c];
                        hi[c] = u;
                        c = p;
                    }
                } else if (v < hd[0]) {
                    hd[0] = v;
                    hi[0] = j;
                    let p = 0;
                    for (;;) {
                        const l = 2 * p + 1;
                        const r = l + 1;
                        let m = p;
                        if (l < n_candidates && hd[l] > hd[m]) m = l;
                        if (r < n_candidates && hd[r] > hd[m]) m = r;
                        if (m === p) break;
                        const t = hd[p];
                        hd[p] = hd[m];
                        hd[m] = t;
                        const u = hi[p];
                        hi[p] = hi[m];
                        hi[m] = u;
                        p = m;
                    }
                }
            }
            if (size < n_candidates) {
                throw new Error(
                    `Only ${size} of the ${n_candidates} requested neighbors exist for element ${i}. ` +
                        `Lower "n_neighbors".`,
                );
            }
            const order = Array.from({ length: size }, (_, c) => c).sort((a, b) => hd[a] - hd[b]);
            for (let c = 0; c < size; ++c) {
                dist[i * n_candidates + c] = hd[order[c]];
                idx[i * n_candidates + c] = hi[order[c]];
            }
        }
        return { dist, idx };
    }

    /**
     * Selects the NN pairs by rescaled distance.
     *
     * Each candidate distance is divided by the local scale of both endpoints — the mean distance
     * to their 4th to 6th neighbors — before the top `n_neighbors` are taken. Without this the NN
     * set is an ordinary KNN graph and dense regions dominate the attractive term.
     *
     * @protected
     * @param {{ dist: Float64Array; idx: Int32Array }} candidates
     * @param {number} n_candidates
     * @param {number} n_neighbors
     * @returns {Int32Array} Flat `[i, j, ...]` pairs.
     */
    _scaled_neighbors(candidates, n_candidates, n_neighbors) {
        const N = this._N;
        const { dist, idx } = candidates;

        // sig is floored: coincident points otherwise give a zero scale and NaN distances.
        const sig = new Float64Array(N);
        for (let i = 0; i < N; ++i) {
            const base = i * n_candidates;
            let sum = 0;
            let count = 0;
            for (let c = 3; c < 6 && c < n_candidates; ++c) {
                sum += dist[base + c];
                ++count;
            }
            sig[i] = Math.max(count > 0 ? sum / count : 0, 1e-10);
        }

        const pairs = new Int32Array(N * n_neighbors * 2);
        const scaled = new Float64Array(n_candidates);
        const order = new Int32Array(n_candidates);
        let p = 0;

        for (let i = 0; i < N; ++i) {
            const base = i * n_candidates;
            for (let c = 0; c < n_candidates; ++c) {
                const j = idx[base + c];
                const dv = dist[base + c];
                scaled[c] = (dv * dv) / (sig[i] * sig[j]);
                order[c] = c;
            }
            order.sort((a, b) => scaled[a] - scaled[b]);
            for (let c = 0; c < n_neighbors; ++c) {
                pairs[p++] = i;
                pairs[p++] = idx[base + order[c]];
            }
        }
        return pairs;
    }

    /**
     * Draws `n_samples` distinct indices, excluding `self` and everything in `reject`.
     *
     * @protected
     * @param {number} n_samples
     * @param {number} self
     * @param {Int32Array | number[]} reject
     * @param {Int32Array} out - Receives the sample; must hold `n_samples`.
     * @returns {number} How many were drawn, short of `n_samples` only if the pool is exhausted.
     */
    _sample_excluding(n_samples, self, reject, out) {
        const N = this._N;
        const randomizer = this._randomizer;
        let count = 0;
        for (let s = 0; s < n_samples; ++s) {
            // The reference retries without limit, which hangs if the pool cannot supply enough
            // distinct points. Bounded here and reported short instead; callers use what they get.
            let attempts = 0;
            for (; attempts < 1000; ++attempts) {
                const j = Math.floor(randomizer.random * N);
                if (j === self || j >= N) continue;
                let taken = false;
                for (let c = 0; c < count && !taken; ++c) if (out[c] === j) taken = true;
                for (let c = 0; c < reject.length && !taken; ++c) if (reject[c] === j) taken = true;
                if (taken) continue;
                out[count++] = j;
                break;
            }
            if (attempts >= 1000) break;
        }
        return count;
    }

    /**
     * Samples mid-near pairs: for each point, six random draws from the whole dataset, of which the
     * second closest is kept. Sampling globally rather than from the neighbour list is what makes
     * these pairs "mid-near" and is the mechanism behind PaCMAP's global structure in phase 1.
     *
     * @protected
     * @param {Matrix} X
     * @param {number} n_MN
     * @returns {Int32Array} Flat `[i, j, ...]` pairs.
     */
    _sample_mn_pairs(X, n_MN) {
        const N = this._N;
        const metric = /** @type {Metric} */ (this.parameter("metric"));
        const pairs = new Int32Array(N * n_MN * 2);
        const picked = new Int32Array(n_MN);
        const sampled = new Int32Array(6);
        const dists = new Float64Array(6);
        let p = 0;

        for (let i = 0; i < N; ++i) {
            const x_i = X.row(i);
            let n_picked = 0;
            for (let m = 0; m < n_MN; ++m) {
                const reject = n_picked > 0 ? picked.subarray(0, n_picked) : EMPTY;
                const got = this._sample_excluding(6, i, reject, sampled);
                if (got < 2) continue;
                for (let c = 0; c < got; ++c) dists[c] = metric(x_i, X.row(sampled[c]));
                // Second closest of the draw: drop the nearest, take the nearest of the rest.
                let best = 0;
                for (let c = 1; c < got; ++c) if (dists[c] < dists[best]) best = c;
                let second = -1;
                for (let c = 0; c < got; ++c) {
                    if (c === best) continue;
                    if (second < 0 || dists[c] < dists[second]) second = c;
                }
                picked[n_picked++] = sampled[second];
                pairs[p++] = i;
                pairs[p++] = sampled[second];
            }
        }
        return pairs.slice(0, p);
    }

    /**
     * Samples further pairs: random points that are not among `i`'s neighbors.
     *
     * @protected
     * @param {Int32Array} nn_pairs
     * @param {number} n_neighbors
     * @param {number} n_FP
     * @returns {Int32Array} Flat `[i, j, ...]` pairs.
     */
    _sample_fp_pairs(nn_pairs, n_neighbors, n_FP) {
        const N = this._N;
        const pairs = new Int32Array(N * n_FP * 2);
        const reject = new Int32Array(n_neighbors);
        const drawn = new Int32Array(n_FP);
        let p = 0;

        for (let i = 0; i < N; ++i) {
            for (let c = 0; c < n_neighbors; ++c) reject[c] = nn_pairs[(i * n_neighbors + c) * 2 + 1];
            const got = this._sample_excluding(n_FP, i, reject, drawn);
            for (let c = 0; c < got; ++c) {
                pairs[p++] = i;
                pairs[p++] = drawn[c];
            }
        }
        return pairs.slice(0, p);
    }

    /**
     * Accumulates the gradient of one pair type into `grad_flat`.
     *
     * All three losses share the form `w * f(d_ij)` with `d_ij = 1 + ||y_i - y_j||²`, so one loop
     * covers them: the attractive pairs use `d_ij / (a + d_ij)` and the repulsive ones
     * `1 / (1 + d_ij)`, which is the same denominator with `a = 1` and the sign flipped.
     *
     * `Y` is indexed flat rather than through `Matrix.row`, which would allocate a subarray per
     * endpoint — two per pair, tens of millions per run, and by far the dominant cost of a step.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N ⨯ d gradient accumulator, modified in place.
     * @param {Int32Array} pairs - Flat `[i, j, ...]` pair array.
     * @param {number} w - Weight for this pair type.
     * @param {number} a - Denominator constant: 10 for NN, 10000 for MN, 1 for FP.
     * @param {boolean} repulsive
     */
    _accumulate_gradients(grad_flat, pairs, w, a, repulsive) {
        if (w === 0) return;
        const Yv = this.Y.values;
        const d = /** @type {number} */ (this.parameter("d"));
        const n_pairs = pairs.length / 2;
        const sign = repulsive ? -1 : 1;

        for (let p = 0; p < n_pairs; ++p) {
            const base_i = pairs[p * 2] * d;
            const base_j = pairs[p * 2 + 1] * d;

            let sq_dist = 0;
            for (let k = 0; k < d; ++k) {
                const diff = Yv[base_i + k] - Yv[base_j + k];
                sq_dist += diff * diff;
            }
            const denom = a + 1 + sq_dist;
            const coeff = (sign * w * 2 * a) / (denom * denom);

            for (let k = 0; k < d; ++k) {
                const g = coeff * (Yv[base_i + k] - Yv[base_j + k]);
                grad_flat[base_i + k] += g;
                grad_flat[base_j + k] -= g;
            }
        }
    }

    /**
     * Returns the weight schedule for the current iteration.
     *
     * @protected
     * @param {number} iter - Current iteration (0-indexed)
     * @returns {{ w_nn: number; w_mn: number; w_fp: number }}
     */
    _get_weights(iter) {
        const num_iters = /** @type {number[]} */ (this.parameter("num_iters"));
        const [p1, p2] = num_iters;
        if (iter < p1) {
            // Phase 1: MN weight linearly decays from 1000 to 3
            const t = iter / p1;
            return { w_nn: 2.0, w_mn: 1000.0 * (1 - t) + 3.0 * t, w_fp: 1.0 };
        } else if (iter < p1 + p2) {
            // Phase 2: fixed weights
            return { w_nn: 3.0, w_mn: 3.0, w_fp: 1.0 };
        } else {
            // Phase 3: MN disabled
            return { w_nn: 1.0, w_mn: 0.0, w_fp: 1.0 };
        }
    }

    /**
     * Applies Adam optimizer update to Y using accumulated gradients.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N ⨯ d gradient
     */
    _adam_update(grad_flat) {
        const lr = /** @type {number} */ (this.parameter("lr"));
        const total = this._N * /** @type {number} */ (this.parameter("d"));
        const beta1 = 0.9;
        const beta2 = 0.999;
        this._adam_t = (this._adam_t ?? 0) + 1;
        const t = /** @type {number} */ (this._adam_t);
        // Bias correction folded into the step size, as in the reference, so `eps` keeps the same
        // weight relative to `sqrt(v)` on every iteration.
        const lr_t = (lr * Math.sqrt(1 - beta2 ** t)) / (1 - beta1 ** t);
        const Yv = this.Y.values;
        const m = /** @type {Float64Array} */ (this._adam_m);
        const v = /** @type {Float64Array} */ (this._adam_v);

        for (let i = 0; i < total; ++i) {
            const g = grad_flat[i];
            m[i] += (1 - beta1) * (g - m[i]);
            v[i] += (1 - beta2) * (g * g - v[i]);
            Yv[i] -= (lr_t * m[i]) / (Math.sqrt(v[i]) + 1e-7);
        }
    }

    /**
     * Initializes PaCMAP: preprocessing, PCA embedding, NN, MN and FP pairs, and Adam state.
     *
     * @returns {this}
     */
    init() {
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const seed = /** @type {number} */ (this.parameter("seed"));
        const n_neighbors = /** @type {number} */ (this.parameter("n_neighbors"));
        const MN_ratio = /** @type {number} */ (this.parameter("MN_ratio"));
        const FP_ratio = /** @type {number} */ (this.parameter("FP_ratio"));
        const apply_pca = /** @type {boolean} */ (this.parameter("apply_pca"));

        // 1. Preprocessing. Above 100 dimensions the reference centers and reduces to 100; below it
        // rescales to the unit range and centers. Neither is cosmetic: the loss constants 10, 10000
        // and 1 are absolute, so an unnormalized input silently changes what "near" means.
        let X_knn;
        let pca_solution = false;
        if (this._D > 100 && apply_pca) {
            X_knn = /** @type {Matrix} */ (PCA.transform(this.X, { d: 100, seed }));
            pca_solution = true;
        } else {
            X_knn = this.X.clone();
            const v = X_knn.values;
            let min = Infinity;
            let max = -Infinity;
            for (let i = 0; i < v.length; ++i) {
                if (v[i] < min) min = v[i];
                if (v[i] > max) max = v[i];
            }
            const range = max - min;
            const scale = range > 0 ? 1 / range : 1;
            for (let i = 0; i < v.length; ++i) v[i] = (v[i] - min) * scale;
            const D = this._D;
            for (let j = 0; j < D; ++j) {
                let mean = 0;
                for (let i = 0; i < N; ++i) mean += v[i * D + j];
                mean /= N;
                for (let i = 0; i < N; ++i) v[i * D + j] -= mean;
            }
        }
        this._X_knn = X_knn;

        // 2. PCA initialization scaled by 0.01.
        const init_source = pca_solution ? X_knn : /** @type {Matrix} */ (PCA.transform(X_knn, { d, seed }));
        this.Y = new Matrix(N, d, (i, j) => init_source.entry(i, j) * 0.01);

        // 3. Pair construction. The 50 extra candidates are what the density rescaling re-ranks.
        const n_candidates = Math.min(n_neighbors + 50, N - 1);
        const candidates = this._find_candidates(X_knn, n_candidates);
        this._nn_pairs = this._scaled_neighbors(candidates, n_candidates, n_neighbors);

        const n_MN = Math.max(1, Math.round(n_neighbors * MN_ratio));
        const n_FP = Math.max(1, Math.round(n_neighbors * FP_ratio));
        this._mn_pairs = this._sample_mn_pairs(X_knn, n_MN);
        this._fp_pairs = this._sample_fp_pairs(this._nn_pairs, n_neighbors, n_FP);

        // 4. Adam optimizer state.
        this._adam_m = new Float64Array(N * d);
        this._adam_v = new Float64Array(N * d);
        this._adam_t = 0;
        this._grad = new Float64Array(N * d);

        this._iter = 0;
        this._is_initialized = true;
        return this;
    }

    /**
     * Accumulates the nearest-neighbor gradient for the current step.
     *
     * Split out from {@link next} because it is the only part of a step {@link LocalMAP} changes,
     * and only in phase 3. Everything else — the mid-near and further terms, the Adam update, the
     * iteration bookkeeping — is shared.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N ⨯ d gradient accumulator, modified in place.
     * @param {number} w_nn
     */
    _accumulate_nn_gradients(grad_flat, w_nn) {
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._nn_pairs), w_nn, 10, false);
    }

    /**
     * Hook run after the embedding has been updated, with the iteration that just finished.
     *
     * Nothing to do here; {@link LocalMAP} redraws its further pairs from it.
     *
     * @protected
     * @param {number} iter - The 0-indexed iteration that just completed.
     */
    _after_step(iter) {
        void iter;
    }

    /**
     * Performs one optimization step.
     *
     * @returns {Matrix}
     */
    next() {
        if (!this._nn_pairs) throw new Error("Call init() first!");
        const iter = this._iter;
        const { w_nn, w_mn, w_fp } = this._get_weights(iter);

        const grad_flat = /** @type {Float64Array} */ (this._grad);
        grad_flat.fill(0);
        this._accumulate_nn_gradients(grad_flat, w_nn);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._mn_pairs), w_mn, 10000, false);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._fp_pairs), w_fp, 1, true);
        this._adam_update(grad_flat);

        this._iter++;
        this._after_step(iter);
        return this.Y;
    }

    /**
     * @param {number} [iterations] - Total number of iterations. Defaults to sum of `num_iters`.
     * @returns {T}
     */
    transform(iterations) {
        const num_iters = /** @type {number[]} */ (this.parameter("num_iters"));
        const total = iterations ?? num_iters.reduce((a, b) => a + b, 0);
        this.check_init();
        for (let i = 0; i < total; ++i) {
            this.next();
        }
        return this.projection;
    }

    /**
     * @param {number} [iterations] - Total number of iterations. Defaults to sum of `num_iters`.
     * @returns {Generator<T, T, void>}
     */
    *generator(iterations) {
        const num_iters = /** @type {number[]} */ (this.parameter("num_iters"));
        const total = iterations ?? num_iters.reduce((a, b) => a + b, 0);
        this.check_init();
        for (let i = 0; i < total; ++i) {
            this.next();
            yield this.projection;
        }
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPaCMAP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new PaCMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPaCMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new PaCMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPaCMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new PaCMAP(X, parameters);
        return dr.transform_async();
    }
}
