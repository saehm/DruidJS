import { euclidean_squared } from "../metrics/euclidean_squared.js";

/**
 * Computes the euclidean distance (`l_2`) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The euclidean distance between `a` and `b`.
 */
export function euclidean(a, b) {
    return Math.sqrt(euclidean_squared(a, b));
}
