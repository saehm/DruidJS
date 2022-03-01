/**
 * @class
 * @alias PCA
 * @augments DR
 */
export class PCA extends DR {
    static principal_components(X: any, parameters: any): Matrix;
    /**
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias PCA
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @returns {PCA}
     */
    constructor(X: Matrix | Array<Array<number>>, parameters: {
        d?: number;
        seed?: number;
        eig_args?: number;
    });
    Y: any[] | Matrix;
    /**
     * Computes the {@link d} principal components of Matrix {@link X}.
     * @returns {Matrix}
     */
    principal_components(): Matrix;
    V: Matrix;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=PCA.d.ts.map