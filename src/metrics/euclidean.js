import { euclidean_squared } from "../metrics/index.js";
/**
 * Computes the euclidean distance (<code>l<sub>2</sub></code>) between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias euclidean
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the euclidean distance between <code>a</code> and <code>b</code>.
 */
export default function (a, b) {
    return Math.sqrt(euclidean_squared(a, b));
}
