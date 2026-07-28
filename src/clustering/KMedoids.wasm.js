/**
 * WASM wrappers for the `KMedoids` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `KMedoids.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * Assigns points to nearest medoids in k-Medoids.
 *
 * Deliberately *not* on a persistent buffer session, despite copying an N ⨯ N distance matrix to do
 * only O(n·k) work. `init()` runs the PAM loop through `_iteration()` and calls this once at the
 * end, so there is no repeat copy to save — only the async `generator()` calls it per iteration. A
 * session here would retain an N ⨯ N buffer for the life of the process and buy nothing.
 *
 * @param {Float64Array} D_val
 * @param {Int32Array} medoids
 * @param {Int32Array} assignments
 * @param {number} n
 * @param {number} k
 * @returns {boolean}
 */
export function wasmKMedoidsAssign(D_val, medoids, assignments, n, k) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrD = alloc(exports, ptrs, D_val.byteLength);
        const ptrM = alloc(exports, ptrs, medoids.byteLength);
        const ptrA = alloc(exports, ptrs, assignments.byteLength);

        new Float64Array(memory.buffer, ptrD, D_val.length).set(D_val);
        new Int32Array(memory.buffer, ptrM, medoids.length).set(medoids);

        exports.kmedoids_assign_f64(ptrD, ptrM, ptrA, n, k);

        assignments.set(new Int32Array(memory.buffer, ptrA, assignments.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}
