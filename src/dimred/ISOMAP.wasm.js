/**
 * WASM wrappers for the `ISOMAP` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `ISOMAP.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * Computes All-Pairs Shortest Path geodesic distances for ISOMAP using WASM Min-Heap Repeated Dijkstra.
 *
 * @param {Int32Array} neighbors_val
 * @param {Float64Array} distances_val
 * @param {Float64Array} out_d_val
 * @param {number} n
 * @param {number} k
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmDijkstraAPSP(neighbors_val, distances_val, out_d_val, n, k) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrN = alloc(exports, ptrs, neighbors_val.byteLength);
        const ptrDist = alloc(exports, ptrs, distances_val.byteLength);
        const ptrOut = alloc(exports, ptrs, out_d_val.byteLength);

        new Int32Array(memory.buffer, ptrN, neighbors_val.length).set(neighbors_val);
        new Float64Array(memory.buffer, ptrDist, distances_val.length).set(distances_val);

        exports.dijkstra_apsp_f64(ptrN, ptrDist, ptrOut, n, k);

        out_d_val.set(new Float64Array(memory.buffer, ptrOut, out_d_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Range-based All-Pairs Shortest Path calculation for worker thread execution.
 *
 * @param {Int32Array} neighbors_val
 * @param {Float64Array} distances_val
 * @param {Float64Array} out_d_val
 * @param {number} n
 * @param {number} k
 * @param {number} start_src
 * @param {number} end_src
 * @returns {boolean}
 */
export function wasmDijkstraAPSPRange(neighbors_val, distances_val, out_d_val, n, k, start_src, end_src) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrN = alloc(exports, ptrs, neighbors_val.byteLength);
        const ptrDist = alloc(exports, ptrs, distances_val.byteLength);
        const ptrOut = alloc(exports, ptrs, out_d_val.byteLength);

        new Int32Array(memory.buffer, ptrN, neighbors_val.length).set(neighbors_val);
        new Float64Array(memory.buffer, ptrDist, distances_val.length).set(distances_val);

        exports.dijkstra_apsp_range_f64(ptrN, ptrDist, ptrOut, n, k, start_src, end_src);

        // Only rows `[start_src, end_src)` were computed. Copying the whole buffer back would
        // overwrite the rows other workers wrote into the same shared output with this instance's
        // untouched (zero) memory.
        const offset = start_src * n;
        const count = (end_src - start_src) * n;
        out_d_val.set(new Float64Array(memory.buffer, ptrOut + offset * 8, count), offset);
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}
