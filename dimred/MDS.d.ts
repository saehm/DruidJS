/**
 * @class
 * @alias MDS
 * @extends DR
 */
export class MDS extends DR {
    /**
     * Classical MDS.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias MDS
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     */
    constructor(X: Matrix, parameters: {
        d?: number;
        metric?: Function | "precomputed";
        seed?: number;
        eig_args?: number;
    });
    _d_X: Matrix;
    Y: Matrix;
    /**
     * @returns {Number} - the stress of the projection.
     */
    stress(): number;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=MDS.d.ts.map