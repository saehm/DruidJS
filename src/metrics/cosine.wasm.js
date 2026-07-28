/**
 * WASM wrappers for the `cosine` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `cosine.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * Computes cosine distance between two 1D vectors A and B.
 *
 * @param {Float64Array | number[]} A
 * @param {Float64Array | number[]} B
 * @returns {number | null} Distance result, or null if WASM is unavailable.
 */
export function wasmCosineDistance(A, B) {
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

        const memA = new Float64Array(memory.buffer, ptrA, len);
        const memB = new Float64Array(memory.buffer, ptrB, len);

        memA.set(A);
        memB.set(B);

        return exports.cosine_distance_f64(ptrA, ptrB, len);
    } finally {
        free_all(exports, ptrs);
    }
}
