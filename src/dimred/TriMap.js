import { spatial_tree } from "../knn/index.js";
import { linspace, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";
import { PCA } from "./PCA.js";
import { wasmTriMapGrad, wasmTriMapUpdate } from "./TriMap.wasm.js";

/** @import {InputType} from "../index.js" */
/** @import {Metric} from "../metrics/index.js" */
/** @import {ParametersTriMap} from "./index.js" */
/** @import {KNN} from "../knn/KNN.js" */

/**
 * The embedding starts as a heavily shrunk PCA projection. Every squared distance in the loss is
 * offset by 1, so an unscaled projection would swamp that term and flatten the gradients.
 */
const INIT_PCA_SCALE = 0.01;
/** Random triplets carry a tenth of the weight of the neighbour-derived ones. */
const RAND_WEIGHT_SCALE = 0.1;
/** Neighbours ranked beyond `n_inliers` that are still considered when sampling. */
const NUM_EXTRA_KNN = 50;
/** Iteration after which the momentum switches from `INIT_MOMENTUM` to `FINAL_MOMENTUM`. */
const SWITCH_ITER_MOMENTUM = 250;
const INIT_MOMENTUM = 0.5;
const FINAL_MOMENTUM = 0.8;

/**
 * Tempered logarithm with temperature `t`, which degrades to `Math.log` as `t` approaches 1.
 *
 * @param {number} x
 * @param {number} t - Temperature.
 * @returns {number}
 */
function tempered_log(x, t) {
    if (Math.abs(t - 1) < 1e-5) return Math.log(x);
    return (1 / (1 - t)) * (x ** (1 - t) - 1);
}

/**
 * TriMap
 *
 * A dimensionality reduction technique that preserves both local and global
 * structure using triplets. It is designed to be a more robust alternative
 * to t-SNE and UMAP.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersTriMap>
 * @category Dimensionality Reduction
 */
export class TriMap extends DR {
    /** @protected */
    get _wasm_session_keys() {
        return ["trimap_grad", "trimap_update"];
    }

    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersTriMap>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
     * @see {@link https://github.com/eamid/trimap}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                weight_temp: 0.5,
                n_inliers: 12,
                n_outliers: 4,
                n_random: 3,
                d: 2,
                metric: euclidean,
                lr: 0.1,
                seed: 1212,
            },
            parameters,
        );
    }

    /**
     * @param {Matrix | null} [pca=null] - Initial Embedding (if null then PCA gets used). Default is `null`
     * @param {import("../knn/KNN.js").KNN<number[] | Float64Array, any> | null} [knn=null] - KNN Object (if null then a KDTree or BallTree gets used, depending on the metric). Default is `null`
     */
    init(pca = null, knn = null) {
        const X = this.X;
        const N = X.shape[0];
        //const c = /** @type {number} */ (this._parameters.c);
        const d = /** @type {number} */ (this._parameters.d);
        const metric = /** @type {Metric} */ (this._parameters.metric);
        const seed = /** @type {number} */ (this._parameters.seed);
        this.n_inliers = /** @type {number} */ (this._parameters.n_inliers);
        this.n_outliers = /** @type {number} */ (this._parameters.n_outliers);
        this.n_random = /** @type {number} */ (this._parameters.n_random);
        this.Y = pca ?? /** @type {Matrix} */ (PCA.transform(X, { d, seed })).mult(INIT_PCA_SCALE);
        this.knn = knn ?? spatial_tree(X.to2dArray(), { metric, seed });
        const { triplets, weights } = this._generate_triplets(this.n_inliers, this.n_outliers, this.n_random);
        this.triplets = triplets;
        this.weights = weights;
        // Triplets are point indices held in a Float64 Matrix; the WASM kernel wants them as i32.
        // They never change after init, so convert once instead of on every gradient step.
        /** @type {Int32Array} */
        this._triplets_int32 = Int32Array.from(triplets.values);
        this.lr = /** @type {number} */ (this.parameter("lr"));
        this.C = Infinity;
        this.vel = new Matrix(N, d, 0);
        this.gain = new Matrix(N, d, 1);
        return this;
    }

    /**
     * Generates {@link n_inliers} x {@link n_outliers} x {@link n_random} triplets.
     *
     * @param {number} n_inliers
     * @param {number} n_outliers
     * @param {number} n_random
     */
    _generate_triplets(n_inliers, n_outliers, n_random) {
        const metric = /** @type {Metric} */ (this._parameters.metric);
        const weight_temp = /** @type {number} */ (this._parameters.weight_temp);
        const X = this.X;
        const N = X.shape[0];
        const knn = this.knn;
        if (!knn) throw new Error("Call init() first!");
        const n_extra = Math.min(n_inliers + NUM_EXTRA_KNN, N);
        const nbrs = new Matrix(N, n_extra);
        const knn_distances = new Matrix(N, n_extra);
        for (let i = 0; i < N; ++i) {
            // Excludes `i` itself and comes back closest first. Filtering on a zero distance instead
            // would also drop genuine duplicates of `i`, which shifts every column that `sig` reads.
            knn.search_by_index(i, n_extra).forEach((d, j) => {
                nbrs.set_entry(i, j, d.index);
                knn_distances.set_entry(i, j, d.distance);
            });
        }
        // Scale parameter: the mean distance to the 3rd, 4th and 5th nearest neighbour. The
        // reference reads columns 3..5 of a row that begins with the point itself at distance 0;
        // these rows exclude it, so the same three neighbours sit one column earlier.
        const sig = new Float64Array(N);
        for (let i = 0; i < N; ++i) {
            sig[i] = Math.max(
                (knn_distances.entry(i, 2) + knn_distances.entry(i, 3) + knn_distances.entry(i, 4)) / 3,
                1e-10,
            );
        }

        const P = this._find_p(knn_distances, sig, nbrs);

        let triplets = this._sample_knn_triplets(P, nbrs, n_inliers, n_outliers);
        let n_triplets = triplets.shape[0];
        const outlier_distances = new Float64Array(n_triplets);
        for (let i = 0; i < n_triplets; ++i) {
            const j = triplets.entry(i, 0);
            const k = triplets.entry(i, 2);
            outlier_distances[i] = metric(X.row(j), X.row(k));
        }
        let weights = this._find_weights(triplets, P, nbrs, outlier_distances, sig);

        if (n_random > 0) {
            const { random_triplets, random_weights } = this._sample_random_triplets(X, n_random, sig);
            triplets = triplets.concat(random_triplets, "vertical");
            weights = Float64Array.from([...weights, ...random_weights.map((w) => RAND_WEIGHT_SCALE * w)]);
        }
        n_triplets = triplets.shape[0];

        // Weights are log-ratios, so they are shifted to be non-negative rather than divided by
        // their maximum, and then compressed with a tempered log.
        let min_weight = Infinity;
        for (let i = 0; i < n_triplets; ++i) {
            if (Number.isNaN(weights[i])) weights[i] = 0;
            if (weights[i] < min_weight) min_weight = weights[i];
        }
        for (let i = 0; i < n_triplets; ++i) {
            weights[i] = tempered_log(1 + (weights[i] - min_weight), weight_temp);
        }
        return {
            triplets: triplets,
            weights: weights,
        };
    }

    /**
     * Calculates the log-similarity matrix P.
     *
     * Kept in log space: exponentiating here underflows to zero for anything but the closest
     * neighbours, and every consumer either compares entries or takes their difference, both of
     * which are order-preserving under the logarithm.
     *
     * @private
     * @param {Matrix} knn_distances - Matrix of pairwise knn distances
     * @param {Float64Array} sig - Scaling factor for the distances
     * @param {Matrix} nbrs - Nearest neighbors
     * @returns {Matrix} Pairwise log-similarity matrix
     */
    _find_p(knn_distances, sig, nbrs) {
        const [N, n_neighbors] = knn_distances.shape;
        return new Matrix(N, n_neighbors, (i, j) => {
            return -(knn_distances.entry(i, j) ** 2) / (sig[i] * sig[nbrs.entry(i, j)]);
        });
    }

    /**
     * Sample nearest neighbors triplets based on the similarity values given in P.
     *
     * @private
     * @param {Matrix} P - Matrix of pairwise similarities between each point and its neighbors given in matrix nbrs.
     * @param {Matrix} nbrs - Nearest neighbors indices for each point. The similarity values are given in matrix
     *   {@link P}. Row i corresponds to the i-th point.
     * @param {number} n_inliers - Number of inlier points.
     * @param {number} n_outliers - Number of outlier points.
     */
    _sample_knn_triplets(P, nbrs, n_inliers, n_outliers) {
        const N = nbrs.shape[0];
        const triplets_list = [];
        for (let i = 0; i < N; ++i) {
            const sort_indices = this.__argsort(P.row(i));
            for (let j = 0; j < n_inliers; ++j) {
                // `sort_indices` holds column offsets into `nbrs`, which already excludes `i`, so
                // the j-th most similar neighbour is taken directly.
                const sim = nbrs.entry(i, sort_indices[j]);
                const rejects = [i, ...Array.from(sort_indices.slice(0, j + 2)).map((idx) => nbrs.entry(i, idx))];
                const samples = this._rejection_sample(n_outliers, N, rejects);
                for (let k = 0; k < samples.length; ++k) {
                    const out = samples[k];
                    triplets_list.push([i, sim, out]);
                }
            }
        }
        const triplets = new Matrix(triplets_list.length, 3);
        for (let t = 0; t < triplets_list.length; ++t) {
            triplets.set_entry(t, 0, triplets_list[t][0]);
            triplets.set_entry(t, 1, triplets_list[t][1]);
            triplets.set_entry(t, 2, triplets_list[t][2]);
        }
        return triplets;
    }

    /**
     * Should do the same as np.argsort()
     *
     * @private
     * @param {Float64Array | number[]} A
     */
    __argsort(A) {
        return linspace(0, A.length - 1).sort((i, j) => A[j] - A[i]);
    }

    /**
     * Samples {@link n_samples} integers from a given interval [0, {@link max_int}] while rejection the values that are
     * in the {@link rejects}.
     *
     * @private
     * @param {number} n_samples
     * @param {number} max_int
     * @param {number[]} rejects
     */
    _rejection_sample(n_samples, max_int, rejects) {
        const randomizer = this._randomizer;
        const interval = linspace(0, max_int - 1).filter((d) => rejects.indexOf(d) < 0);
        return randomizer.choice(interval, Math.min(n_samples, interval.length));
    }

    /**
     * Calculates the weights for the sampled nearest neighbors triplets
     *
     * @private
     * @param {Matrix} triplets - Sampled Triplets.
     * @param {Matrix} P - Pairwise similarity matrix.
     * @param {Matrix} nbrs - Nearest Neighbors
     * @param {Float64Array} outlier_distances - Matrix of pairwise outlier distances
     * @param {Float64Array} sig - Scaling factor for the distances.
     */
    _find_weights(triplets, P, nbrs, outlier_distances, sig) {
        const n_triplets = triplets.shape[0];
        const weights = new Float64Array(n_triplets);
        for (let t = 0; t < n_triplets; ++t) {
            const i = triplets.entry(t, 0);
            const sim = nbrs.row(i).indexOf(triplets.entry(t, 1));
            const p_sim = P.entry(i, sim);
            const p_out = -(outlier_distances[t] ** 2) / (sig[i] * sig[triplets.entry(t, 2)]);
            // Both are log-similarities, so their ratio is the difference.
            weights[t] = p_sim - p_out;
        }
        return weights;
    }

    /**
     * Sample uniformly ranom triplets
     *
     * @private
     * @param {Matrix} X - Data matrix.
     * @param {number} n_random - Number of random triplets per point
     * @param {Float64Array} sig - Scaling factor for the distances
     */
    _sample_random_triplets(X, n_random, sig) {
        const metric = /** @type {Metric} */ (this.parameter("metric"));
        const randomizer = this._randomizer;
        const N = X.shape[0];
        const random_triplets = new Matrix(N * n_random, 3);
        const random_weights = new Float64Array(N * n_random);
        for (let i = 0; i < N; ++i) {
            const n_i = i * n_random;
            const indices = Array.from({ length: N }, (_, idx) => idx).filter((idx) => idx !== i);
            for (let j = 0; j < n_random; ++j) {
                let [sim, out] = randomizer.choice(indices, 2);
                let p_sim = -(metric(X.row(i), X.row(sim)) ** 2) / (sig[i] * sig[sim]);
                let p_out = -(metric(X.row(i), X.row(out)) ** 2) / (sig[i] * sig[out]);

                // The nearer of the pair has to play the inlier; comparing the logs orders them
                // the same way comparing the similarities would.
                if (p_sim < p_out) {
                    [sim, out] = [out, sim];
                    [p_sim, p_out] = [p_out, p_sim];
                }
                const index = n_i + j;
                random_triplets.set_entry(index, 0, i);
                random_triplets.set_entry(index, 1, sim);
                random_triplets.set_entry(index, 2, out);
                random_weights[index] = p_sim - p_out;
            }
        }
        return {
            random_triplets: random_triplets,
            random_weights: random_weights,
        };
    }

    /**
     * Computes the gradient for updating the embedding.
     *
     * @param {Matrix} Y - The embedding
     * @returns {{ grad: Matrix; loss: number }} The gradient and the current loss.
     */
    _grad(Y) {
        const n_inliers = this.n_inliers;
        const n_outliers = this.n_outliers;
        const triplets = this.triplets;
        const weights = this.weights;
        if (!triplets || n_inliers === undefined || n_outliers === undefined || !weights || !this._triplets_int32)
            throw new Error("Call init() first!");
        const [N, dim] = Y.shape;

        const grad = new Matrix(N, dim, 0);
        const loss = wasmTriMapGrad(
            Y.values,
            this._triplets_int32,
            weights,
            grad.values,
            N,
            dim,
            n_inliers,
            n_outliers,
        );

        if (loss !== null) {
            return { grad, loss };
        }

        const n_triplets = triplets.shape[0];
        const y_ij = new Float64Array(dim);
        const y_ik = new Float64Array(dim);
        let d_ij = 1;
        let d_ik = 1;
        let js_loss = 0;
        const n_knn_triplets = N * n_inliers * n_outliers;

        for (let t = 0; t < n_triplets; ++t) {
            const [i, j, k] = triplets.row(t);
            if (t % n_outliers === 0 || t >= n_knn_triplets) {
                d_ij = 1;
                d_ik = 1;
                for (let d = 0; d < dim; ++d) {
                    const Y_id = Y.entry(i, d);
                    const Y_jd = Y.entry(j, d);
                    const Y_kd = Y.entry(k, d);
                    y_ij[d] = Y_id - Y_jd;
                    y_ik[d] = Y_id - Y_kd;
                    d_ij += y_ij[d] ** 2;
                    d_ik += y_ik[d] ** 2;
                }
            } else {
                d_ik = 1;
                for (let d = 0; d < dim; ++d) {
                    const Y_id = Y.entry(i, d);
                    const Y_kd = Y.entry(k, d);
                    y_ik[d] = Y_id - Y_kd;
                    d_ik += y_ik[d] ** 2;
                }
            }

            js_loss += weights[t] / (1 + d_ik / d_ij);
            const w = weights[t] / (d_ij + d_ik) ** 2;
            for (let d = 0; d < dim; ++d) {
                const gs = y_ij[d] * d_ik * w;
                const go = y_ik[d] * d_ij * w;
                grad.add_entry(i, d, gs - go);
                grad.sub_entry(j, d, gs);
                grad.add_entry(k, d, go);
            }
        }
        return { grad, loss: js_loss };
    }

    /**
     * @param {number} max_iteration
     * @returns {T}
     */
    transform(max_iteration = 400) {
        this.check_init();
        try {
            for (let iter = 0; iter < max_iteration; ++iter) {
                this._next(iter);
            }
            return this.projection;
        } finally {
            this._release_wasm();
        }
    }

    /**
     * @param {number} max_iteration
     * @returns {Generator<T, T, void>}
     */
    *generator(max_iteration = 400) {
        this.check_init();
        try {
            for (let iter = 0; iter < max_iteration; ++iter) {
                this._next(iter);
                yield this.projection;
            }
            return this.projection;
        } finally {
            this._release_wasm();
        }
    }

    /**
     * Does the iteration step.
     *
     * @private
     * @param {number} iter
     */
    _next(iter) {
        const gamma = iter > SWITCH_ITER_MOMENTUM ? FINAL_MOMENTUM : INIT_MOMENTUM;
        const vel = this.vel;
        const triplets = this.triplets;
        if (!vel || !triplets || this.lr === undefined) throw new Error("Call init() first!");
        const Y = this.Y.add(vel.mult(gamma));
        const { grad, loss } = this._grad(Y);
        // Reported per triplet, so the value does not scale with the triplet count.
        this.C = loss / triplets.shape[0];
        this.Y = this._update_embedding(Y, iter, grad);
        // The delta-bar-delta optimizer adapts a per-coordinate gain, so the global learning rate
        // is left alone; only the plain gradient-descent variants anneal it.
        return this.Y;
    }

    /**
     * Updates the embedding.
     *
     * @private
     * @param {Matrix} Y
     * @param {number} iter
     * @param {Matrix} grad
     */
    _update_embedding(Y, iter, grad) {
        const [N, dim] = Y.shape;
        const gamma = iter > SWITCH_ITER_MOMENTUM ? FINAL_MOMENTUM : INIT_MOMENTUM;
        const gain = this.gain;
        const vel = this.vel;
        const lr = this.lr;
        if (!vel || !gain || lr === undefined) throw new Error("Call init() first!");

        if (wasmTriMapUpdate(Y.values, grad.values, vel.values, gain.values, N, dim, gamma, lr)) {
            return Y;
        }

        const min_gain = 0.01;
        for (let i = 0; i < N; ++i) {
            for (let d = 0; d < dim; ++d) {
                const new_gain =
                    Math.sign(vel.entry(i, d)) !== Math.sign(grad.entry(i, d))
                        ? gain.entry(i, d) + 0.2
                        : Math.max(gain.entry(i, d) * 0.8, min_gain);
                gain.set_entry(i, d, new_gain);
                vel.set_entry(i, d, gamma * vel.entry(i, d) - lr * gain.entry(i, d) * grad.entry(i, d));
                Y.set_entry(i, d, Y.entry(i, d) + vel.entry(i, d));
            }
        }
        return Y;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTriMap>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new TriMap(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTriMap>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new TriMap(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersTriMap>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new TriMap(X, parameters);
        return dr.transform_async();
    }
}
