/**
 * @class
 * @alias ISOMAP
 * @extends DR
 */
export class ISOMAP extends DR {
    /**
     * Isometric feature mapping (ISOMAP).
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias ISOMAP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} parameters.neighbors - the number of neighbors {@link ISOMAP} should use to project the data.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://doi.org/10.1126/science.290.5500.2319}
     */
    constructor(X: Matrix, parameters: {
        neighbors: number;
        d?: number;
        metric?: Function;
        seed?: number;
        eig_args?: number;
    });
    Y: Matrix;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=ISOMAP.d.ts.map