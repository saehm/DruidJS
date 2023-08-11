import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";

/**
 * @class
 * @alias FASTMAP
 * @extends DR
 */
export class FASTMAP extends DR {
    /**
     * FastMap: a fast algorithm for indexing, data-mining and visualization of traditional and multimedia datasets
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias FASTMAP
     * @param {Matrix} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the dimensionality of the projection.
     * @returns {FASTMAP}
     * @see {@link https://doi.org/10.1145/223784.223812}
     */
    constructor(X, parameters) {
        super(X, { d: 2, metric: euclidean, seed: 1212 }, parameters);
        return this;
    }

    /**
     * Chooses two points which are the most distant in the actual projection.
     * @private
     * @param {function} dist
     * @returns {number[]} An array consisting of first index, second index, and distance between the two points.
     */
    _choose_distant_objects(dist) {
        const X = this.X;
        const N = X.shape[0];
        let a_index = (this._randomizer.random_int % N) - 1;
        let b_index = null;
        let max_dist = -Infinity;
        for (let i = 0; i < N; ++i) {
            const d_ai = dist(a_index, i);
            if (d_ai > max_dist) {
                max_dist = d_ai;
                b_index = i;
            }
        }
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
     * @returns {Matrix} The {@link d}-dimensional projection of the data matrix {@link X}.
     */
    transform() {
        const X = this.X;
        const N = X.shape[0];
        const { d, metric } = this._parameters;
        const Y = new Matrix(N, d, 0);
        let dist = (a, b) => metric(X.row(a), X.row(b));

        for (let _col = 0; _col < d; ++_col) {
            let old_dist = dist;
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
}
