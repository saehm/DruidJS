import { WASM_MIN_SIMPLE_VECTOR_LENGTH } from "../wasm/thresholds.js";
import { wasmChebyshevDistance } from "./chebyshev.wasm.js";

/**
 * Computes the chebyshev distance (L<sub>∞</sub>) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The chebyshev distance between `a` and `b`.
 */
export function chebyshev(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;

    if (n >= WASM_MIN_SIMPLE_VECTOR_LENGTH) {
        const wasmRes = wasmChebyshevDistance(a, b);
        if (wasmRes !== null) {
            return wasmRes;
        }
    }

    const res = [];
    for (let i = 0; i < n; ++i) {
        res.push(Math.abs(a[i] - b[i]));
    }
    return Math.max(...res);
}
