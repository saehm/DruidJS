/**
 * Computes the cosine distance (not similarity) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The cosine distance between `a` and `b`.
 * @example
 * import { cosine } from "@saehrimnir/druidjs";
 * const a = [1, 2, 3];
 * const b = [4, 5, 6];
 * const distance = cosine(a, b); // 0.9746318461970762
 */
export function cosine(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
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
