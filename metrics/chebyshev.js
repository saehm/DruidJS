/**
 * Computes the chebyshev distance (L<sub>∞</sub>) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias chebyshev
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the chebyshev distance between {@link a} and {@link b}.
 */
export default function (a, b) {
    if (a.length != b.length) return undefined;
    let n = a.length;
    let res = [];
    for (let i = 0; i < n; ++i) {
        res.push(Math.abs(a[i] - b[i]));
    }
    return Math.max(...res);
}
