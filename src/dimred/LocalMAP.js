import { PaCMAP } from "./PaCMAP.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersLocalMAP} from "./index.js" */

/**
 * LocalMAP
 *
 * A variant of PaCMAP that improves local cluster separation in phase 3 by:
 * 1. Applying a local scaling factor (nn_scale / sqrt(d_ij)) to NN gradients,
 *    amplifying attraction for already-close NN pairs and dampening it for far ones.
 * 2. Periodically resampling FP pairs from random candidates within a distance
 *    threshold in the current embedding, rather than using static random pairs.
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
     * In phase 3, applies local NN scaling and resamples FP pairs every 10 iterations
     * from random candidates within the distance threshold.
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

        // Phase 3: local NN scaling + periodic FP resampling within distance threshold
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const { w_nn, w_mn, w_fp } = this._get_weights(this._iter);
        const nn_scale = (this._low_dist_thres ?? 10) / 2;

        // Resample FP pairs from random embedding-local candidates every 10 iterations
        // (skip the very first phase 3 step to match reference behaviour)
        if (this._iter > phase3_start && this._iter % 10 === 0) {
            this._resample_local_fp_pairs();
        }

        const grad_flat = new Float64Array(N * d);
        // Local NN: attractive with nn_scale/sqrt(d_ij) distance-based scaling
        this._accumulate_gradients_local_nn(grad_flat, /** @type {Int32Array} */ (this._nn_pairs), w_nn, nn_scale);
        // MN: no-op in phase 3 (w_mn = 0), kept for correctness
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._mn_pairs), w_mn, 10000, false);
        // FP: standard repulsive on (periodically resampled) local pairs
        this._accumulate_gradients(grad_flat, /** @type {Int32Array} */ (this._fp_pairs), w_fp, 0, true);

        this._adam_update(grad_flat);
        this._iter++;
        return this.Y;
    }

    /**
     * Accumulates NN gradients with LocalMAP's local scaling: nn_scale / sqrt(d_ij).
     * Amplifies attraction for NN pairs already close in embedding; dampens it for far pairs.
     *
     * @private
     * @param {Float64Array} grad_flat - Flat N×d gradient accumulator (modified in place)
     * @param {Int32Array} pairs - Flat [i, j, i, j, ...] pair array
     * @param {number} w_nn - NN weight
     * @param {number} nn_scale - low_dist_thres / 2
     */
    _accumulate_gradients_local_nn(grad_flat, pairs, w_nn, nn_scale) {
        if (w_nn === 0) return;
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

            // NN loss gradient scaled by nn_scale / sqrt(d_ij)
            const denom = 10 + d_ij;
            const coeff = (w_nn * 20 * nn_scale) / (denom * denom * Math.sqrt(d_ij));

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
     * Resamples FP pairs in place using random candidates within the embedding distance
     * threshold. For each pair (i, j), tries up to 64 random points; if one falls within
     * low_dist_thres and is not a high-dim neighbor of i, it replaces j. Otherwise j is kept.
     *
     * @private
     */
    _resample_local_fp_pairs() {
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const low_dist_thres = this._low_dist_thres ?? 10;
        const threshold_sq = low_dist_thres * low_dist_thres;
        const nn_sets = this._nn_sets_cache;
        const fp_pairs = /** @type {Int32Array} */ (this._fp_pairs);
        const n_pairs = fp_pairs.length / 2;
        const Y = this.Y;
        const randomizer = this._randomizer;

        for (let p = 0; p < n_pairs; ++p) {
            const i = fp_pairs[p * 2];
            const nn_set = nn_sets ? nn_sets[i] : null;
            const y_i = Y.row(i);

            for (let s = 0; s < 64; ++s) {
                const j = randomizer.random_int % N;
                if (j === i || (nn_set && nn_set.has(j))) continue;

                const y_j = Y.row(j);
                let sq_dist = 0;
                for (let k = 0; k < d; ++k) {
                    const diff = y_i[k] - y_j[k];
                    sq_dist += diff * diff;
                }
                if (sq_dist <= threshold_sq) {
                    fp_pairs[p * 2 + 1] = j;
                    break;
                }
            }
            // If no valid candidate found, keep the existing pair
        }
    }

    /**
     * Initializes LocalMAP (same as PaCMAP, but caches nn_sets for phase 3 FP resampling).
     *
     * @returns {LocalMAP<T>}
     */
    init() {
        super.init();
        // Cache low_dist_thres from sealed parameters (avoids type indexing issues)
        this._low_dist_thres = /** @type {number} */ (/** @type {any} */ (this._parameters)["low_dist_thres"] ?? 10);
        // Cache nn_sets for use in phase 3 FP resampling
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
