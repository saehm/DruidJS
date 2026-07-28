/**
 * WASM wrappers for the `KMeans` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `KMeans.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * Assigns points to nearest centroids in k-Means.
 *
 * @param {Float64Array} X_val
 * @param {Float64Array} C_val
 * @param {Int32Array} assignments
 * @param {number} n
 * @param {number} k
 * @param {number} d
 * @returns {boolean} True if executed successfully via WASM.
 */
export function wasmKMeansAssign(X_val, C_val, assignments, n, k, d) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);
        const ptrA = alloc(exports, ptrs, assignments.byteLength);

        const memX = new Float64Array(memory.buffer, ptrX, X_val.length);
        const memC = new Float64Array(memory.buffer, ptrC, C_val.length);
        memX.set(X_val);
        memC.set(C_val);

        exports.kmeans_assign_f64(ptrX, ptrC, ptrA, n, k, d);

        const memA = new Int32Array(memory.buffer, ptrA, assignments.length);
        assignments.set(memA);
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}
