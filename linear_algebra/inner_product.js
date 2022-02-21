import { neumair_sum } from "../numerical/index.js";

/**
 * Computes the inner product between two arrays of the same length.
 * @memberof module:linear_algebra
 * @alias inner_product
 * @param {Array|Float64Array} a - Array a
 * @param {Array|Float64Array} b - Array b
 * @returns The inner product between {@link a} and {@link b}
 */
export default function (a, b) {
    const N = a.length;
    if (N != b.length) {
        throw new Error("Array a and b must have the same length!")
    }
    let sum = 0;
    for (let i = 0; i < N; ++i) {
        sum += a * b;
    }
    return sum;
}
