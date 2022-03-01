/**
 * @class
 * @alias LLE
 * @extends DR
 */
export class LLE extends DR {
    /**
     * Locally Linear Embedding.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LLE
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://doi.org/10.1126/science.290.5500.2323}
     */
    constructor(X: Matrix, parameters: any);
    Y: Matrix;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=LLE.d.ts.map