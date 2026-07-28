/**
 * WASM wrappers for the `Matrix` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `Matrix.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * Computes C = A * B using WASM SIMD.
 * A is (rows_A x cols_A), B is (cols_A x cols_B), out_data is (rows_A x cols_B).
 *
 * @param {Float64Array} A_val
 * @param {number} rows_A
 * @param {number} cols_A
 * @param {Float64Array} B_val
 * @param {number} cols_B
 * @param {Float64Array} C_val
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmMatMul(A_val, rows_A, cols_A, B_val, cols_B, C_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    const sizeA = A_val.byteLength;
    const sizeB = B_val.byteLength;
    const sizeC = C_val.byteLength;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, sizeA);
        const ptrB = alloc(exports, ptrs, sizeB);
        const ptrC = alloc(exports, ptrs, sizeC);

        const memA = new Float64Array(memory.buffer, ptrA, A_val.length);
        memA.set(A_val);

        const memB = new Float64Array(memory.buffer, ptrB, B_val.length);
        memB.set(B_val);

        exports.matmul_f64(ptrA, ptrB, ptrC, rows_A, cols_A, cols_B);

        C_val.set(new Float64Array(memory.buffer, ptrC, C_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Range-based Parallel Matrix Multiplication for worker thread execution.
 *
 * @param {Float64Array} A_val
 * @param {number} rows_A
 * @param {number} cols_A
 * @param {Float64Array} B_val
 * @param {number} cols_B
 * @param {Float64Array} C_val
 * @param {number} start_row
 * @param {number} end_row
 * @returns {boolean}
 */
export function wasmMatMulRange(A_val, rows_A, cols_A, B_val, cols_B, C_val, start_row, end_row) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrB = alloc(exports, ptrs, B_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrB, B_val.length).set(B_val);

        exports.matmul_range_f64(ptrA, ptrB, ptrC, rows_A, cols_A, cols_B, start_row, end_row);

        // Slice only, so concurrent workers sharing `C_val` do not clobber each other's rows.
        const offset = start_row * cols_B;
        const count = (end_row - start_row) * cols_B;
        C_val.set(new Float64Array(memory.buffer, ptrC + offset * 8, count), offset);
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes C = A^T * B using WASM SIMD.
 * A is (cols_A x rows_A), B is (cols_A x cols_B), out_data is (rows_A x cols_B).
 *
 * @param {Float64Array} A_val
 * @param {number} cols_A
 * @param {number} rows_A
 * @param {Float64Array} B_val
 * @param {number} cols_B
 * @param {Float64Array} C_val
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmTransDot(A_val, cols_A, rows_A, B_val, cols_B, C_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrB = alloc(exports, ptrs, B_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrB, B_val.length).set(B_val);

        exports.transDot_f64(ptrA, ptrB, ptrC, cols_A, rows_A, cols_B);

        C_val.set(new Float64Array(memory.buffer, ptrC, C_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes C = A * B^T using WASM SIMD.
 * A is (rows_A x cols_A), B is (cols_B x cols_A), out_data is (rows_A x cols_B).
 *
 * @param {Float64Array} A_val
 * @param {number} rows_A
 * @param {number} cols_A
 * @param {Float64Array} B_val
 * @param {number} cols_B
 * @param {Float64Array} C_val
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmDotTrans(A_val, rows_A, cols_A, B_val, cols_B, C_val) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrB = alloc(exports, ptrs, B_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrB, B_val.length).set(B_val);

        exports.dotTrans_f64(ptrA, ptrB, ptrC, rows_A, cols_A, cols_B);

        C_val.set(new Float64Array(memory.buffer, ptrC, C_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes outer product C = A x B using WASM SIMD.
 *
 * @param {Float64Array} A_val
 * @param {Float64Array} B_val
 * @param {Float64Array} C_val
 * @param {number} len
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmOuter(A_val, B_val, C_val, len) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrB = alloc(exports, ptrs, B_val.byteLength);
        const ptrC = alloc(exports, ptrs, C_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrB, B_val.length).set(B_val);

        exports.outer_f64(ptrA, ptrB, ptrC, len);

        C_val.set(new Float64Array(memory.buffer, ptrC, C_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Computes Y = A * X using WASM SIMD.
 *
 * @param {Float64Array} A_val
 * @param {Float64Array} X_val
 * @param {Float64Array} Y_val
 * @param {number} rows
 * @param {number} cols
 * @returns {boolean}
 */
export function wasmMatVecMul(A_val, X_val, Y_val, rows, cols) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrA = alloc(exports, ptrs, A_val.byteLength);
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);

        new Float64Array(memory.buffer, ptrA, A_val.length).set(A_val);
        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);

        exports.mat_vec_mul_f64(ptrA, ptrX, ptrY, rows, cols);

        Y_val.set(new Float64Array(memory.buffer, ptrY, Y_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}
