/**
 * @class
 * @alias LSP
 * @extends DR
 */
export class LSP extends DR {
    /**
     * Least Squares Projection.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LSP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.neighbors = Math.max(Math.floor(N / 10), 2)] - number of neighbors to consider.
     * @param {Number} [parameters.control_points = Math.ceil(Math.sqrt(N))] - number of controlpoints
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {LSP}
     * @see {@link https://ieeexplore.ieee.org/document/4378370}
     * @todo accept precomputed distance matrix.
     */
    constructor(X: Matrix, parameters: {
        neighbors?: number;
        control_points?: number;
        d?: number;
        metric?: Function;
        seed?: number;
    });
    /**
     *
     * @param {DR} DR - method used for position control points.
     * @param {Object} DR_parameters - Object containing parameters for the DR method which projects the control points
     * @returns {LSP}
     */
    init(DR?: DR, DR_parameters?: any, KNN?: typeof BallTree): LSP;
    _A: Matrix;
    _b: Matrix;
    Y: Matrix;
}
import { DR } from "./DR.js";
import { BallTree } from "../knn/index.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=LSP.d.ts.map