/**
 * @class
 * @alias SAMMON
 * @extends DR
 */
export class SAMMON extends DR {
    /**
     * SAMMON's Mapping
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias SAMMON
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {"PCA"|"MDS"|"random"} [parameters.init = "random"] - Either "PCA" or "MDS", with which SAMMON initialiates the projection. With "random" a random matrix gets used as starting point.
     * @param {Object} [parameters.init_parameters] - Parameters for the {@link init}-DR method.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {SAMMON}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
    constructor(X: Matrix, parameters: {
        d?: number;
        metric?: Function | "precomputed";
        init?: "PCA" | "MDS" | "random";
        init_parameters?: any;
        seed?: number;
    });
    /**
     * initializes the projection.
     * @private
     */
    private init;
    Y: Matrix;
    distance_matrix: Matrix;
    _step(): Matrix;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=SAMMON.d.ts.map