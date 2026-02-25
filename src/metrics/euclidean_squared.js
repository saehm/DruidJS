/**
 * Computes the squared euclidean distance (l_2^2) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The squared euclidean distance between `a` and `b`.

 */
export function euclidean_squared(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        const a_b = a[i] - b[i];
        sum += a_b * a_b;
    }
    return sum;
}
