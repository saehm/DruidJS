import { BallTree } from "../knn/index.js";
import { Matrix } from "../matrix/index.js";
import { euclidean, euclidean_squared } from "../metrics/index.js";
import { DR } from "./DR.js";
import { PCA } from "./PCA.js";

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {ParametersPaCMAP} from "./index.js" */

/**
 * Pairwise Controlled Manifold Approximation Projection (PaCMAP)
 *
 * A dimensionality reduction technique that uses three types of point pairs —
 * nearest neighbor (NN), mid-near (MN), and further (FP) pairs — with a
 * dynamic three-phase weight schedule and Adam optimization to preserve both
 * local and global structure.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersPaCMAP>
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
     * @param {Partial<ParametersPaCMAP>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters) {
        super(
            X,
            {
                n_neighbors: 10,
                MN_ratio: 0.5,
                FP_ratio: 2.0,
                d: 2,
                metric: euclidean,
                lr: 1.0,
                num_iters: [100, 100, 250],
                seed: 1212,
            },
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
     * Samples mid-near pairs for each point.
     * For each point i, repeats n_MN times: samples 6 random non-neighbor
     * candidates, picks the 2nd closest by high-dim distance.
     *
     * @protected
     * @param {Set<number>[]} nn_sets - Array of neighbor index sets per point
     * @param {number} n_MN - Number of mid-near pairs per point
     * @returns {Int32Array} Flat array of [i, j] pairs
     */
    _sample_mn_pairs(nn_sets, n_MN) {
        const N = this._N;
        const X = this.X;
        const randomizer = this._randomizer;
        const pairs = new Int32Array(N * n_MN * 2);
        let idx = 0;

        for (let i = 0; i < N; ++i) {
            const nn_set = nn_sets[i];
            const x_i = X.row(i);

            for (let k = 0; k < n_MN; ++k) {
                // Sample 6 random non-neighbor candidates
                /** @type {{idx: number, dist: number}[]} */
                const candidates = [];
                let attempts = 0;
                while (candidates.length < 6 && attempts < N * 2) {
                    const j = randomizer.random_int % N;
                    attempts++;
                    if (j !== i && !nn_set.has(j)) {
                        candidates.push({ idx: j, dist: euclidean_squared(x_i, X.row(j)) });
                    }
                }
                if (candidates.length < 2) {
                    // Fallback: use any available non-self point
                    const j = (i + 1) % N;
                    pairs[idx++] = i;
                    pairs[idx++] = j;
                    continue;
                }
                candidates.sort((a, b) => a.dist - b.dist);
                // Pick the 2nd closest (index 1)
                pairs[idx++] = i;
                pairs[idx++] = candidates[1].idx;
            }
        }
        return pairs.slice(0, idx);
    }

    /**
     * Samples further pairs for each point (random non-neighbors).
     *
     * @protected
     * @param {Set<number>[]} nn_sets - Array of neighbor index sets per point
     * @param {number} n_FP - Number of further pairs per point
     * @returns {Int32Array} Flat array of [i, j] pairs
     */
    _sample_fp_pairs(nn_sets, n_FP) {
        const N = this._N;
        const randomizer = this._randomizer;
        const pairs = new Int32Array(N * n_FP * 2);
        let idx = 0;

        for (let i = 0; i < N; ++i) {
            const nn_set = nn_sets[i];
            let count = 0;
            let attempts = 0;
            while (count < n_FP && attempts < N * 3) {
                const j = randomizer.random_int % N;
                attempts++;
                if (j !== i && !nn_set.has(j)) {
                    pairs[idx++] = i;
                    pairs[idx++] = j;
                    count++;
                }
            }
        }
        return pairs.slice(0, idx);
    }

    /**
     * Computes gradient coefficients and updates the gradient matrix for one pair type.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N×d gradient accumulator (modified in place)
     * @param {Int32Array} pairs - Flat [i, j, i, j, ...] pair array
     * @param {number} w - Weight for this pair type
     * @param {number} attr_num - Numerator constant for attractive (10 for NN, 10000 for MN); 0 for repulsive
     * @param {boolean} repulsive - Whether this is a repulsive pair type
     */
    _accumulate_gradients(grad_flat, pairs, w, attr_num, repulsive) {
        if (w === 0) return;
        const Y = this.Y;
        const d = /** @type {number} */ (this.parameter("d"));
        const n_pairs = pairs.length / 2;

        for (let p = 0; p < n_pairs; ++p) {
            const i = pairs[p * 2];
            const j = pairs[p * 2 + 1];
            const y_i = Y.row(i);
            const y_j = Y.row(j);

            // d_ij = 1 + ||y_i - y_j||²
            let sq_dist = 0;
            for (let k = 0; k < d; ++k) {
                const diff = y_i[k] - y_j[k];
                sq_dist += diff * diff;
            }
            const d_ij = 1 + sq_dist;

            /** @type {number} */
            let coeff;
            if (repulsive) {
                // FP loss: 1/(1+d_ij), gradient: -2/(1+d_ij)²
                coeff = (-w * 2) / (d_ij * d_ij);
            } else {
                // NN loss: d_ij/(attr_num+d_ij), gradient: 2*attr_num/(attr_num+d_ij)²
                const denom = attr_num + d_ij;
                coeff = (w * 2 * attr_num) / (denom * denom);
            }

            const base_i = i * d;
            const base_j = j * d;
            for (let k = 0; k < d; ++k) {
                const diff = y_i[k] - y_j[k];
                const g = coeff * diff;
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
     * @param {Float64Array} grad_flat - Flat N×d gradient
     */
    _adam_update(grad_flat) {
        const lr = /** @type {number} */ (this.parameter("lr"));
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const beta1 = 0.9;
        const beta2 = 0.999;
        const eps = 1e-7;
        this._adam_t = (this._adam_t ?? 0) + 1;
        const t = /** @type {number} */ (this._adam_t);
        const bc1 = 1 - beta1 ** t;
        const bc2 = 1 - beta2 ** t;
        const Y = this.Y;
        const m = /** @type {Float64Array} */ (this._adam_m);
        const v = /** @type {Float64Array} */ (this._adam_v);

        for (let i = 0; i < N; ++i) {
            const base = i * d;
            const y_i = Y.row(i);
            for (let k = 0; k < d; ++k) {
                const g = grad_flat[base + k];
                m[base + k] = beta1 * m[base + k] + (1 - beta1) * g;
                v[base + k] = beta2 * v[base + k] + (1 - beta2) * g * g;
                const m_hat = m[base + k] / bc1;
                const v_hat = v[base + k] / bc2;
                y_i[k] -= lr * (m_hat / (Math.sqrt(v_hat) + eps));
            }
        }
    }

    /**
     * Initializes PaCMAP: PCA embedding, KNN pairs, MN pairs, FP pairs, Adam state.
     *
     * @returns {PaCMAP<T>}
     */
    init() {
        const X = this.X;
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const seed = /** @type {number} */ (this.parameter("seed"));
        const metric = /** @type {Metric} */ (this.parameter("metric"));
        const n_neighbors = /** @type {number} */ (this.parameter("n_neighbors"));
        const MN_ratio = /** @type {number} */ (this.parameter("MN_ratio"));
        const FP_ratio = /** @type {number} */ (this.parameter("FP_ratio"));

        // 1. PCA initialization scaled by 0.01 (X is always Matrix here)
        const pca_init = /** @type {Matrix} */ (PCA.transform(X, { d, seed }));
        this.Y = new Matrix(N, d, (i, j) => pca_init.entry(i, j) * 0.01);

        // 2. Build KNN graph for NN pairs
        const knn = new BallTree(X.to2dArray(), { metric, seed });
        const n_MN = Math.max(1, Math.round(n_neighbors * MN_ratio));
        const n_FP = Math.max(1, Math.round(n_neighbors * FP_ratio));

        /** @type {Int32Array[]} */
        const nn_indices = [];
        /** @type {Set<number>[]} */
        const nn_sets = [];
        // NN pairs: flat [i, j] pairs
        const nn_pairs = new Int32Array(N * n_neighbors * 2);
        let nn_idx = 0;

        for (let i = 0; i < N; ++i) {
            const neighbors = knn.search(X.row(i), n_neighbors + 1);
            /** @type {number[]} */
            const idxs = [];
            for (const nb of neighbors) {
                if (nb.index !== i) idxs.push(nb.index);
                if (idxs.length >= n_neighbors) break;
            }
            nn_indices[i] = new Int32Array(idxs);
            nn_sets[i] = new Set(idxs);
            for (const j of idxs) {
                nn_pairs[nn_idx++] = i;
                nn_pairs[nn_idx++] = j;
            }
        }
        this._nn_pairs = nn_pairs.slice(0, nn_idx);

        // 3. MN pairs (mid-near sampling)
        this._mn_pairs = this._sample_mn_pairs(nn_sets, n_MN);

        // 4. FP pairs (random non-neighbors)
        this._fp_pairs = this._sample_fp_pairs(nn_sets, n_FP);

        // 5. Adam optimizer state
        this._adam_m = new Float64Array(N * d);
        this._adam_v = new Float64Array(N * d);
        this._adam_t = 0;

        this._iter = 0;
        return this;
    }

    /**
     * Performs one optimization step.
     *
     * @returns {Matrix}
     */
    next() {
        if (!this._nn_pairs) throw new Error("Call init() first!");
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const { w_nn, w_mn, w_fp } = this._get_weights(this._iter);

        const grad_flat = new Float64Array(N * d);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._nn_pairs), w_nn, 10, false);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._mn_pairs), w_mn, 10000, false);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._fp_pairs), w_fp, 0, true);
        this._adam_update(grad_flat);

        this._iter++;
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
