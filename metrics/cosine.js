/**
 * Computes the cosine distance (not similarity) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias cosine
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} The cosine distance between {@link a} and {@link b}.
 * 
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 * 
 * druid.cosine([1,0],[1,1]) == 0.7853981633974484 == π/4;
 * 
 */
export default function (a, b) {
    if (a.length !== b.length) return undefined;
    let n = a.length;
    let sum = 0;
    let sum_a = 0;
    let sum_b = 0;
    for (let i = 0; i < n; ++i) {
        sum += a[i] * b[i];
        sum_a += a[i] * a[i];
        sum_b += b[i] * b[i];
    }
    return Math.acos(sum / (Math.sqrt(sum_a) * Math.sqrt(sum_b)));
}
