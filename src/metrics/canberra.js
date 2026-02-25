/**
 * Computes the canberra distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The canberra distance between `a` and `b`.
 * @see {@link https://en.wikipedia.org/wiki/Canberra_distance}
 */
export function canberra(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += Math.abs(a[i] - b[i]) / (Math.abs(a[i]) + Math.abs(b[i]));
    }
    return sum;
}
