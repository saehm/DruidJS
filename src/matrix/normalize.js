import { euclidean } from "../metrics/index.js";
import { norm } from "./index.js";

/** @import { Metric } from "../metrics/index.js" */

/**
 * Normalizes Vector `v`.
 *
 * @category Matrix
 * @param {number[] | Float64Array} v - Vector
 * @param {Metric} metric
 * @returns {number[] | Float64Array} - The normalized vector with length 1.
 */
export function normalize(v, metric = euclidean) {
    const v_norm = norm(v, metric);
    return v.map((value) => value / v_norm);
}
