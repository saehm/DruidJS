import { simultaneous_poweriteration} from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";

/**
 * @class
 * @alias HIPP
 * @extends DR
 */
export class HIPP{
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias HIPP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @returns {HIPP}
     */
    constructor(X, d=2) {
        return this;
    }


    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    transform() {
        return Y;
    }

    /**
     * @returns {Matrix} Returns the projection.
     */
    get projection() {
        return this.Y;
    }
} 