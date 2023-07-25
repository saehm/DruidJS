/**
 * Computes the squared euclidean distance (l<sub>2</sub><sup>2</sup>) between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias euclidean_squared
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the squared euclidean distance between <code>a</code> and <code>b</code>.
 */
export default function (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        const a_b = a[i] - b[i];
        sum += a_b * a_b;
    }
    return sum;
}
