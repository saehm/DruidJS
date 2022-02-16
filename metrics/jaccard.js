/**
 * Computes the jaccard distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias jaccard
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the jaccard distance between {@link a} and {@link b}.
 */
export default function (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    let num_non_zero = 0;
    let num_equal = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] != 0;
        const y = b[i] != 0;
        num_non_zero += x || y;
        num_equal += x && y;
    }
    return (num_non_zero - num_equal) / num_non_zero;
}
