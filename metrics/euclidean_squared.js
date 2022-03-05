import { neumair_sum } from "../numerical/index.js";
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
    const s = new Float64Array(n);
    for (let i = 0; i < n; ++i) {
        const x = a[i];
        const y = b[i];
        const x_y = x - y;
        s[i] = x_y * x_y;
    }
    return neumair_sum(s);
}
