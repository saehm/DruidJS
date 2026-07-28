/**
 * WASM wrappers for the `normalize` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `normalize.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * In-place Vector Normalization using WASM SIMD.
 *
 * @param {Float64Array} V_val
 * @returns {boolean} True if executed via WASM.
 */
export function wasmNormalize(V_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = V_val.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrV = alloc(exports, ptrs, V_val.byteLength);

        new Float64Array(memory.buffer, ptrV, len).set(V_val);
        exports.normalize_f64(ptrV, len);
        V_val.set(new Float64Array(memory.buffer, ptrV, len));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}
