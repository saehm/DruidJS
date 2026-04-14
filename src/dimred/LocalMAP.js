import { BallTree } from "../knn/index.js";
import { euclidean_squared } from "../metrics/index.js";
import { PaCMAP } from "./PaCMAP.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersLocalMAP} from "./index.js" */

/**
 * LocalMAP
 *
 * A variant of PaCMAP that improves local cluster separation by dynamically
 * resampling further pairs (FP) in phase 3 using nearby points in the current
 * low-dimensional embedding space, rather than randomly sampled non-neighbors.
 *
 * @class
 * @template {InputType} T
 * @extends PaCMAP<T>
 * @category Dimensionality Reduction
 * @see {@link https://arxiv.org/abs/2012.04456|PaCMAP Paper}
 * @see {@link PaCMAP} for the base algorithm
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const X = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]];
 * const localmap = new druid.LocalMAP(X, {
 *     n_neighbors: 10,
 *     low_dist_thres: 10,
 *     seed: 42
 * });
 *
 * const Y = localmap.transform(); // 450 iterations (default)
 * // [[x1, y1], [x2, y2], [x3, y3]]
 */
export class LocalMAP extends PaCMAP {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLocalMAP>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters = {}) {
        // Merge low_dist_thres into parameters before the seal in DR constructor.
        // DR.constructor does Object.seal({ ...defaults, ...parameters }), so
        // passing low_dist_thres here ensures it lands in the sealed object.
        super(X, { low_dist_thres: 10, ...parameters });
    }

    /**
     * Performs one optimization step.
     * In phase 3, resamples FP pairs from the current embedding space
     * and applies distance-based weight scaling.
     *
     * @returns {import("../matrix/index.js").Matrix}
     */
    next() {
        if (!this._nn_pairs) throw new Error("Call init() first!");
        const num_iters = /** @type {number[]} */ (this.parameter("num_iters"));
        const phase3_start = num_iters[0] + num_iters[1];

        if (this._iter < phase3_start) {
            // Phases 1 and 2: identical to PaCMAP
            return super.next();
        }

        // Phase 3: resample FP pairs from embedding neighbors
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const low_dist_thres = /** @type {number} */ (this._low_dist_thres ?? 10);
        const low_dist_thres_sq = low_dist_thres * low_dist_thres;
        const { w_nn, w_mn, w_fp } = this._get_weights(this._iter);
        const Y = this.Y;

        // Build a KNN structure on the current embedding to find nearby pairs
        const y_array = Y.to2dArray();
        const seed = /** @type {number} */ (this.parameter("seed"));
        const emb_knn = new BallTree(y_array, { metric: euclidean_squared, seed });
        const n_FP = /** @type {number} */ (this.parameter("FP_ratio")) *
            /** @type {number} */ (this.parameter("n_neighbors"));
        const n_FP_int = Math.max(1, Math.round(n_FP));

        // Build local FP pairs by finding nearby embedding neighbors
        // (these are not NN neighbors in high-dim space but nearby in low-dim)
        const nn_sets = this._nn_sets_cache;
        /** @type {Int32Array} */
        let local_fp_pairs;

        if (nn_sets) {
            const pair_buf = new Int32Array(N * n_FP_int * 2);
            let idx = 0;
            for (let i = 0; i < N; ++i) {
                const nn_set = nn_sets[i];
                // Search for nearby points in embedding space
                const neighbors = emb_knn.search(y_array[i], Math.min(n_FP_int * 3 + 1, N));
                let count = 0;
                for (const nb of neighbors) {
                    if (nb.index === i || (nn_set && nn_set.has(nb.index))) continue;
                    pair_buf[idx++] = i;
                    pair_buf[idx++] = nb.index;
                    count++;
                    if (count >= n_FP_int) break;
                }
            }
            local_fp_pairs = pair_buf.slice(0, idx);
        } else {
            local_fp_pairs = /** @type {Int32Array} */ (this._fp_pairs);
        }

        // Accumulate gradients with local weight scaling for FP pairs
        const grad_flat = new Float64Array(N * d);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._nn_pairs), w_nn, 10, false);
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._mn_pairs), w_mn, 10000, false);

        // FP pairs with local distance-based weight scaling
        this._accumulate_gradients_local_fp(
            grad_flat,
            local_fp_pairs,
            w_fp,
            low_dist_thres,
            low_dist_thres_sq,
        );

        this._adam_update(grad_flat);
        this._iter++;
        return this.Y;
    }

    /**
     * Accumulates FP gradients with LocalMAP's distance-based weight scaling.
     * For pairs within low_dist_thres, scales w_fp by low_dist_thres / (2 * sqrt(d_ij)).
     *
     * @private
     * @param {Float64Array} grad_flat - Flat N×d gradient accumulator (modified in place)
     * @param {Int32Array} pairs - Flat [i, j, i, j, ...] pair array
     * @param {number} w_fp - Base FP weight
     * @param {number} low_dist_thres - Distance threshold
     * @param {number} low_dist_thres_sq - Squared distance threshold
     */
    _accumulate_gradients_local_fp(grad_flat, pairs, w_fp, low_dist_thres, low_dist_thres_sq) {
        if (w_fp === 0) return;
        const Y = this.Y;
        const d = /** @type {number} */ (this.parameter("d"));
        const n_pairs = pairs.length / 2;

        for (let p = 0; p < n_pairs; ++p) {
            const i = pairs[p * 2];
            const j = pairs[p * 2 + 1];
            const y_i = Y.row(i);
            const y_j = Y.row(j);

            let sq_dist = 0;
            for (let k = 0; k < d; ++k) {
                const diff = y_i[k] - y_j[k];
                sq_dist += diff * diff;
            }
            const d_ij = 1 + sq_dist;

            // Apply local weight scaling when pair is within distance threshold
            let w = w_fp;
            if (sq_dist < low_dist_thres_sq) {
                w *= low_dist_thres / (2 * Math.sqrt(d_ij));
            }

            // FP loss: 1/(1+d_ij), gradient: -2/(1+d_ij)²
            const coeff = (-w * 2) / (d_ij * d_ij);

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
     * Initializes LocalMAP (same as PaCMAP, but caches nn_sets for phase 3 resampling).
     *
     * @returns {LocalMAP<T>}
     */
    init() {
        super.init();
        // Cache low_dist_thres from sealed parameters (avoids type indexing issues)
        this._low_dist_thres = /** @type {number} */ (/** @type {any} */ (this._parameters)["low_dist_thres"] ?? 10);
        // Cache nn_sets for use in phase 3 FP resampling
        // We rebuild them from _nn_pairs
        const N = this._N;
        const nn_pairs = this._nn_pairs;
        if (!nn_pairs) return this;
        /** @type {Set<number>[]} */
        const nn_sets = Array.from({ length: N }, () => new Set());
        for (let p = 0; p < nn_pairs.length; p += 2) {
            nn_sets[nn_pairs[p]].add(nn_pairs[p + 1]);
        }
        this._nn_sets_cache = nn_sets;
        return this;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLocalMAP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new LocalMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLocalMAP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new LocalMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLocalMAP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new LocalMAP(X, parameters);
        return dr.transform_async();
    }
}
