/**
 * WASM wrappers for the `neumair_sum` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `neumair_sum.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * Computes Neumaier sum using WASM SIMD.
 *
 * @param {Float64Array | number[]} V_val
 * @returns {number | null}
 */
export function wasmNeumaierSum(V_val) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = V_val.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrV = alloc(exports, ptrs, len * 8);

        new Float64Array(memory.buffer, ptrV, len).set(V_val);
        return exports.neumair_sum_f64(ptrV, len);
    } finally {
        free_all(exports, ptrs);
    }
}
