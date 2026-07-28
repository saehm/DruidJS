import { WASM_MIN_SIMPLE_VECTOR_LENGTH } from "../wasm/thresholds.js";
import { wasmManhattanDistance } from "./manhattan.wasm.js";

/**
 * Computes the manhattan distance (`l_1`) between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The manhattan distance between `a` and `b`.
 */
export function manhattan(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;

    if (n >= WASM_MIN_SIMPLE_VECTOR_LENGTH) {
        const wasmRes = wasmManhattanDistance(a, b);
        if (wasmRes !== null) {
            return wasmRes;
        }
    }

    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += Math.abs(a[i] - b[i]);
    }
    return sum;
}
