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
 * Exact k-nearest neighbors for rows `[start_row, end_row)` using WASM SIMD.
 *
 * Selection happens inside the kernel, so the copy back is `k` entries per row rather than `n`, and
 * the caller can block the search without ever materializing an n ⨯ n distance matrix.
 *
 * @param {Float64Array} X_val - Row-major `n` ⨯ `d` input.
 * @param {number} n
 * @param {number} d
 * @param {number} k - Neighbors per row, excluding the row itself.
 * @param {Float64Array} out_val - `n` ⨯ `2k`: `k` ascending distances then their `k` indices.
 *   Rows outside the range are left untouched.
 * @param {number} start_row
 * @param {number} end_row
 * @returns {boolean} False if WASM is unavailable and the caller must fall back to JS.
 */
export function wasmKnnBlock(X_val, n, d, k, out_val, start_row, end_row) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrOut = alloc(exports, ptrs, out_val.byteLength);

        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);

        exports.euclidean_knn_block_f64(ptrX, ptrOut, n, d, k, start_row, end_row);

        // Slice only, so concurrent workers sharing `out_val` do not clobber each other's rows.
        const row_len = 2 * k;
        const offset = start_row * row_len;
        const count = (end_row - start_row) * row_len;
        out_val.set(new Float64Array(memory.buffer, ptrOut + offset * 8, count), offset);
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
