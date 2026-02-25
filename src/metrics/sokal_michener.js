/**
 * Computes the Sokal-Michener distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The Sokal-Michener distance between `a` and `b`.

 */
export function sokal_michener(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let num_not_equal = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] !== 0;
        const y = b[i] !== 0;
        num_not_equal += x !== y ? 1 : 0;
    }
    return (2 * num_not_equal) / (n + num_not_equal);
}
