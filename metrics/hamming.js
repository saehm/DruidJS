/**
 * Computes the hamming distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias hamming
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the hamming distance between {@link a} and {@link b}.
 */
export default function (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    let disagree = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i];
        const y = b[i];
        disagree += x != y;
    }
    return disagree / n;
}
