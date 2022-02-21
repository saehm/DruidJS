import { norm } from "./index.js";
import { euclidean } from "../metrics/index.js";

/**
 * Normalizes Vector {@link v}.
 * @memberof module:matrix
 * @alias normalize
 * @param {Array<Number>|Float64Array} v - Vector
 * @param {Function} metric 
 * @returns {Array<Number>|Float64Array} - The normalized vector with length 1.
 */
export default function(v, metric = euclidean)  {
    const v_norm = norm(v, metric);
    return v.map(value => value / v_norm);
}