/**
 * WASM wrappers for the `bray_curtis` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `bray_curtis.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * Computes Bray-Curtis distance.
 *
 * @param {Float64Array | number[]} A
 * @param {Float64Array | number[]} B
 * @returns {number | null}
 */
export function wasmBrayCurtisDistance(A, B) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;
    const len = A.length;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, len * 8);
        const ptrB = alloc(exports, ptrs, len * 8);

        new Float64Array(memory.buffer, ptrA, len).set(A);
        new Float64Array(memory.buffer, ptrB, len).set(B);
        return exports.bray_curtis_distance_f64(ptrA, ptrB, len);
    } finally {
        free_all(exports, ptrs);
    }
}
