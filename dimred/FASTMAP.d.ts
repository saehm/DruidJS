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
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the dimensionality of the projection.
     * @returns {FASTMAP}
     * @see {@link https://doi.org/10.1145/223784.223812}
     */
    constructor(X: Matrix, parameters: {
        d?: number;
        metric?: Function;
        seed?: number;
    });
    /**
     * Chooses two points which are the most distant in the actual projection.
     * @private
     * @param {Function} dist
     * @returns {Array} An array consisting of first index, second index, and distance between the two points.
     */
    private _choose_distant_objects;
    Y: Matrix;
}
import { DR } from "./DR.js";
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=FASTMAP.d.ts.map