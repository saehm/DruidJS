import { neumair_sum } from "../numerical/index.js";
/**
 * Computes the squared euclidean distance (l<sub>2</sub><sup>2</sup>) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias euclidean_squared
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the squared euclidean distance between {@link a} and {@link b}.
 */
export default function (a, b) {
    if (a.length != b.length) return undefined;
    let n = a.length;
    let s = new Array(n);
    for (let i = 0; i < n; ++i) {
        let x = a[i];
        let y = b[i];
        s[i] = (x - y) * (x - y);
    }
    return neumair_sum(s);
}
