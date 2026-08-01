import { PaCMAP } from "./PaCMAP.js";

/** @import {InputType} from "../index.js" */
/** @import {Matrix} from "../matrix/index.js" */
/** @import {ParametersLocalMAP} from "./index.js" */

/**
 * LocalMAP
 *
 * A variant of {@link PaCMAP} that sharpens local cluster separation in phase 3. Attraction on the
 * nearest-neighbor pairs is scaled by `low_dist_thres / (2 * sqrt(d_ij))`, which strengthens it for
 * pairs already close in the embedding and weakens it for far ones, and the further pairs are
 * periodically redrawn from points that are *near* in the current embedding rather than staying the
 * random set chosen at initialization.
 *
 * @class
 * @template {InputType} T
 * @extends PaCMAP<T, ParametersLocalMAP>
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
        // `DR` seals `{ ...defaults, ...parameters }`, so the extra key has to be present here to
        // survive; it cannot be added to the sealed object afterwards.
        super(X, { low_dist_thres: 10, ...parameters });
    }

    /**
     * The first iteration of phase 3, past which LocalMAP diverges from PaCMAP.
     *
     * The reference tests `itr > phase_1 + phase_2`, so the first phase 3 step still runs plain
     * PaCMAP; this is that boundary, and "after" means strictly greater.
     *
     * @private
     * @returns {number}
     */
    get _phase3_start() {
        const num_iters = /** @type {number[]} */ (this.parameter("num_iters"));
        return num_iters[0] + num_iters[1];
    }

    /**
     * Phases 1 and 2 attract exactly as PaCMAP does; phase 3 applies the local scaling.
     *
     * @protected
     * @param {Float64Array} grad_flat
     * @param {number} w_nn
     */
    _accumulate_nn_gradients(grad_flat, w_nn) {
        if (this._iter <= this._phase3_start) return super._accumulate_nn_gradients(grad_flat, w_nn);
        const low_dist_thres = /** @type {number} */ (this.parameter("low_dist_thres"));
        this._accumulate_gradients_local_nn(
            grad_flat,
            /** @type {Int32Array} */ (this._nn_pairs),
            w_nn,
            low_dist_thres / 2,
        );
    }

    /**
     * Redraws the further pairs every tenth iteration of phase 3.
     *
     * Runs after the embedding update, so the redraw sees the layout the step just produced, and
     * is keyed on the iteration that finished — `itr % 10` in the reference is the loop variable,
     * not the next one.
     *
     * @protected
     * @param {number} iter
     */
    _after_step(iter) {
        if (iter > this._phase3_start && iter % 10 === 0) {
            this._resample_local_fp_pairs(/** @type {number} */ (this.parameter("low_dist_thres")));
        }
    }

    /**
     * Accumulates NN gradients with LocalMAP's local scaling, `nn_scale / sqrt(d_ij)`.
     *
     * @protected
     * @param {Float64Array} grad_flat - Flat N ⨯ d gradient accumulator, modified in place.
     * @param {Int32Array} pairs - Flat `[i, j, ...]` pair array.
     * @param {number} w_nn - NN weight.
     * @param {number} nn_scale - `low_dist_thres / 2`.
     */
    _accumulate_gradients_local_nn(grad_flat, pairs, w_nn, nn_scale) {
        if (w_nn === 0) return;
        const Yv = this.Y.values;
        const d = /** @type {number} */ (this.parameter("d"));
        const n_pairs = pairs.length / 2;

        for (let p = 0; p < n_pairs; ++p) {
            const base_i = pairs[p * 2] * d;
            const base_j = pairs[p * 2 + 1] * d;

            let sq_dist = 0;
            for (let k = 0; k < d; ++k) {
                const diff = Yv[base_i + k] - Yv[base_j + k];
                sq_dist += diff * diff;
            }
            const d_ij = 1 + sq_dist;
            const denom = 10 + d_ij;
            const coeff = (w_nn * 20 * nn_scale) / (denom * denom * Math.sqrt(d_ij));

            for (let k = 0; k < d; ++k) {
                const g = coeff * (Yv[base_i + k] - Yv[base_j + k]);
                grad_flat[base_i + k] += g;
                grad_flat[base_j + k] -= g;
            }
        }
    }

    /**
     * Redraws the further pairs from points within `low_dist_thres` of `i` in the embedding.
     *
     * A point that is far in the input but has drifted close in the embedding is exactly what the
     * repulsive term needs to push apart, so sampling from the current layout rather than the
     * original random set is what sharpens the cluster boundaries. Candidates that are already
     * neighbors are rejected, and a row that finds nothing within the threshold keeps its old
     * partner.
     *
     * @protected
     * @param {number} low_dist_thres
     */
    _resample_local_fp_pairs(low_dist_thres) {
        const N = this._N;
        const d = /** @type {number} */ (this.parameter("d"));
        const n_neighbors = /** @type {number} */ (this.parameter("n_neighbors"));
        const threshold_sq = low_dist_thres * low_dist_thres;
        const nn_pairs = /** @type {Int32Array} */ (this._nn_pairs);
        const fp_pairs = /** @type {Int32Array} */ (this._fp_pairs);
        const n_FP = fp_pairs.length / 2 / N;
        const Yv = this.Y.values;
        const randomizer = this._randomizer;
        const drawn = new Int32Array(n_FP);

        for (let i = 0; i < N; ++i) {
            const base_i = i * d;
            for (let s = 0; s < n_FP; ++s) {
                // -1 marks "gave up", and stays -1 through the rest of the row: the old partner
                // that replaces it is substituted at write-out and so never blocks a later draw,
                // which is what the reference's sentinel does.
                drawn[s] = -1;
                // Deliberate deviation: the reference's give-up counter is only reached by
                // duplicate and neighbor rejections — a candidate rejected for being too far
                // `continue`s past the check — so a point with nothing inside the threshold retries
                // until it happens to draw 100 duplicates, which is unbounded work as N grows.
                // Every draw counts here instead. Measured against the reference across
                // low_dist_thres 1 to 40 on IRIS, the difference is under 3% on cluster separation.
                for (let attempt = 0; attempt < 100; ++attempt) {
                    const j = Math.floor(randomizer.random * N);
                    if (j === i || j >= N) continue;

                    let taken = false;
                    for (let c = 0; c < s && !taken; ++c) if (drawn[c] === j) taken = true;
                    for (let c = 0; c < n_neighbors && !taken; ++c) {
                        if (nn_pairs[(i * n_neighbors + c) * 2 + 1] === j) taken = true;
                    }
                    if (taken) continue;

                    const base_j = j * d;
                    let sq_dist = 0;
                    for (let k = 0; k < d; ++k) {
                        const diff = Yv[base_i + k] - Yv[base_j + k];
                        sq_dist += diff * diff;
                    }
                    if (sq_dist > threshold_sq) continue;

                    drawn[s] = j;
                    break;
                }
            }
            for (let s = 0; s < n_FP; ++s) {
                if (drawn[s] >= 0) fp_pairs[(i * n_FP + s) * 2 + 1] = drawn[s];
            }
        }
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
