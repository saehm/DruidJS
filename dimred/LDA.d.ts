/**
 * @class
 * @alias LDA
 * @extends DR
 */
export class LDA extends DR {
    /**
     * Linear Discriminant Analysis.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LDA
     * @param {Matrix} X - The high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Array} parameters.labels - The labels / classes for each data point.
     * @param {number} [parameters.d = 2] - The dimensionality of the projection.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x}
     */
    constructor(X: Matrix, parameters: {
        labels: any[];
        d?: number;
        seed?: number;
        eig_args?: number;
    });
    Y: any[] | Matrix;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=LDA.d.ts.map