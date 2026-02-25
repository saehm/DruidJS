/**
 * Computes the Bray-Curtis distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The Bray-Curtis distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity}
 */
export function bray_curtis(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    let sum_abs_diff = 0;
    let sum_ab = 0;
    for (let i = 0; i < a.length; ++i) {
        sum_abs_diff += Math.abs(a[i] - b[i]);
        sum_ab += a[i] + b[i];
    }
    return sum_abs_diff / sum_ab;
}
