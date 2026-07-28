/**
 * WASM wrappers for the `distance_matrix` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `distance_matrix.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * Computes Pairwise Euclidean Distance Matrix using WASM SIMD.
 *
 * @param {Float64Array} X_val
 * @param {number} n
 * @param {number} d
 * @param {Float64Array} D_val
 * @returns {boolean}
 */
export function wasmDistanceMatrix(X_val, n, d, D_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrD = alloc(exports, ptrs, D_val.byteLength);

        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);

        exports.euclidean_distance_matrix_f64(ptrX, ptrD, n, d);

        D_val.set(new Float64Array(memory.buffer, ptrD, D_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Range-based Pairwise Euclidean Distance Matrix for worker thread execution.
 *
 * @param {Float64Array} X_val
 * @param {number} n
 * @param {number} d
 * @param {Float64Array} D_val
 * @param {number} start_row
 * @param {number} end_row
 * @returns {boolean}
 */
export function wasmDistanceMatrixRange(X_val, n, d, D_val, start_row, end_row) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrD = alloc(exports, ptrs, D_val.byteLength);

        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);

        exports.euclidean_distance_matrix_range_f64(ptrX, ptrD, n, d, start_row, end_row);

        // Slice only, so concurrent workers sharing `D_val` do not clobber each other's rows.
        const offset = start_row * n;
        const count = (end_row - start_row) * n;
        D_val.set(new Float64Array(memory.buffer, ptrD + offset * 8, count), offset);
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}
