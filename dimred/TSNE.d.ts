/**
 * @class
 * @alias TSNE
 * @extends DR
 */
export class TSNE extends DR {
    /**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TSNE
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.perplexity = 50] - perplexity.
     * @param {Number} [parameters.epsilon = 10] - learning parameter.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TSNE}
     */
    constructor(X: Matrix, parameters: {
        perplexity?: number;
        epsilon?: number;
        d?: number;
        metric?: Function | "precomputed";
        seed?: number;
    });
    _iter: number;
    Y: Matrix;
    /**
     *
     * @param {Matrix} distance_matrix - accepts a precomputed distance matrix
     * @returns {TSNE}
     */
    init(): TSNE;
    _ystep: Matrix;
    _gains: Matrix;
    _P: Matrix;
    /**
     * performs a optimization step
     * @private
     * @returns {Matrix}
     */
    private next;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=TSNE.d.ts.map