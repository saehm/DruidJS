/**
 * @class
 * @alias LTSA
 * @extends DR
 */
export class LTSA extends DR {
    /**
     * Local Tangent Space Alignment
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LTSA
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} parameters.neighbors - the number of neighbors {@link LTSA} should use to project the data.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
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
//# sourceMappingURL=LTSA.d.ts.map