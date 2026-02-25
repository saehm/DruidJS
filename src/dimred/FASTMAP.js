import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";

/** @import { InputType } from "../index.js" */
/** @import { ParametersFASTMAP } from "./index.js"; */

/**
 * FastMap algorithm for dimensionality reduction.
 *
 * A very fast algorithm for projecting high-dimensional data into a lower-dimensional
 * space while preserving pairwise distances. It works similarly to PCA but uses
 * only a subset of the data to find projection axes.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersFASTMAP>
 * @category Dimensionality Reduction
 */
export class FASTMAP extends DR {
    /**
     * FastMap: a fast algorithm for indexing, data-mining and visualization of traditional and multimedia datasets.
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersFASTMAP>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1145/223784.223812}
     */
    constructor(X, parameters) {
        super(X, { d: 2, metric: euclidean, seed: 1212 }, parameters);
    }

    /**
     * Chooses two points which are the most distant in the actual projection.
     *
     * @private
     * @param {(a: number, b: number) => number} dist
     * @returns {[number, number, number]} An array consisting of first index, second index, and distance between the
     *   two points.
     */
    _choose_distant_objects(dist) {
        const X = this.X;
        const N = X.shape[0];
        let a_index = this._randomizer.random_int % N;
        /** @type {number | null} */
        let b_index = null;
        let max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            const d_ai = dist(a_index, i);
            if (d_ai > max_dist) {
                max_dist = d_ai;
                b_index = i;
            }
        }
        if (b_index === null) throw new Error("should not happen!");
        max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            const d_bi = dist(b_index, i);
            if (d_bi > max_dist) {
                max_dist = d_bi;
                a_index = i;
            }
        }
        return [a_index, b_index, max_dist];
    }

    /**
     * Computes the projection.
     *
     * @returns {T} The `d`-dimensional projection of the data matrix `X`.
     */
    transform() {
        const X = this.X;
        const N = X.shape[0];
        const d = /** @type {number} */ (this._parameters.d);
        const metric = /** @type {typeof euclidean} */ (this._parameters.metric);
        const Y = new Matrix(N, d, 0);
        /** @type {(a: number, b: number) => number} */
        let dist = (a, b) => metric(X.row(a), X.row(b));

        for (let _col = 0; _col < d; ++_col) {
            const old_dist = dist;
            // choose pivot objects
            const [a_index, b_index, d_ab] = this._choose_distant_objects(dist);
            if (d_ab !== 0) {
                // project the objects on the line (O_a, O_b)
                for (let i = 0; i < N; ++i) {
                    const d_ai = dist(a_index, i);
                    const d_bi = dist(b_index, i);
                    const y_i = (d_ai ** 2 + d_ab ** 2 - d_bi ** 2) / (2 * d_ab);
                    Y.set_entry(i, _col, y_i);
                }
                // consider the projections of the objects on a
                // hyperplane perpendicluar to the line (a, b);
                // the distance function D'() between two
                // projections is given by Eq.4
                dist = (a, b) => Math.sqrt(old_dist(a, b) ** 2 - (Y.entry(a, _col) - Y.entry(b, _col)) ** 2);
            }
        }
        // return embedding.
        this.Y = Y;
        return this.projection;
    }

    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersFASTMAP>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new FASTMAP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersFASTMAP>} parameters
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new FASTMAP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersFASTMAP>} parameters
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new FASTMAP(X, parameters);
        return dr.transform_async();
    }
}
