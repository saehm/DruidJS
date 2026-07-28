import { WASM_MIN_SIMPLE_VECTOR_LENGTH } from "../wasm/thresholds.js";
import { wasmInnerProduct } from "./inner_product.wasm.js";

/**
 * Computes the inner product between two arrays of the same length.
 *
 * @category Linear Algebra
 * @param {number[] | Float64Array} a - Array a.
 * @param {number[] | Float64Array} b - Array b.
 * @returns The inner product between `a` and `b`.
 */
export function inner_product(a, b) {
    const N = a.length;
    if (N !== b.length) {
        throw new Error("Array a and b must have the same length!");
    }
    if (N >= WASM_MIN_SIMPLE_VECTOR_LENGTH) {
        const wasmRes = wasmInnerProduct(a, b);
        if (wasmRes !== null) {
            return wasmRes;
        }
    }
    let sum = 0;
    for (let i = 0; i < N; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}
