import { simultaneous_poweriteration} from "../linear_algebra/index";
import { Matrix } from "../matrix/index";

/**
 * @class
 * @alias DR
 */
export class DR{
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias DR
     * @param {Matrix} X - the high-dimensional data. 
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @returns {DR}
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