// src/matrix/Matrix.as.ts

import { axpy_simd_f64, dot_simd_f64, zero_f64 } from "../wasm/shared.as";

/**
 * Matrix Multiplication: C (rows_A x cols_B) = A (rows_A x cols_A) * B (cols_A x cols_B)
 */
export function matmul_f64(
    a_ptr: usize,
    b_ptr: usize,
    out_ptr: usize,
    rows_A: i32,
    cols_A: i32,
    cols_B: i32
): void {
    matmul_range_f64(a_ptr, b_ptr, out_ptr, rows_A, cols_A, cols_B, 0, rows_A);
}

/**
 * Range-based Parallel Matrix Multiplication Kernel.
 * Allows worker threads to compute matrix rows [start_row, end_row) in parallel.
 */
export function matmul_range_f64(
    a_ptr: usize,
    b_ptr: usize,
    out_ptr: usize,
    rows_A: i32,
    cols_A: i32,
    cols_B: i32,
    start_row: i32,
    end_row: i32
): void {
    for (let i = start_row; i < end_row; ++i) {
        const i_cols_A = i * cols_A;
        const out_row = out_ptr + ((i * cols_B) << 3);

        zero_f64(out_row, cols_B);

        for (let k = 0; k < cols_A; ++k) {
            const aik = load<f64>(a_ptr + ((i_cols_A + k) << 3));
            if (aik == 0.0) continue;
            axpy_simd_f64(out_row, b_ptr + ((k * cols_B) << 3), aik, cols_B);
        }
    }
}

/**
 * Transposed Matrix Multiplication: C (rows_A x cols_B) = A^T * B
 * A is stored as (cols_A x rows_A) in row-major order.
 *
 * The k loop is outermost, so unlike {@link matmul_range_f64} the output cannot be zeroed per row
 * as it is reached -- every row is touched on every k.
 */
export function transDot_f64(
    a_ptr: usize,
    b_ptr: usize,
    out_ptr: usize,
    cols_A: i32,
    rows_A: i32,
    cols_B: i32
): void {
    zero_f64(out_ptr, rows_A * cols_B);

    for (let k = 0; k < cols_A; ++k) {
        const k_rows_A = k * rows_A;
        const b_row = b_ptr + ((k * cols_B) << 3);
        for (let i = 0; i < rows_A; ++i) {
            const aki = load<f64>(a_ptr + ((k_rows_A + i) << 3));
            if (aki == 0.0) continue;
            axpy_simd_f64(out_ptr + ((i * cols_B) << 3), b_row, aki, cols_B);
        }
    }
}

/**
 * Dot product with transposed B: C (rows_A x cols_B) = A * B^T
 * B is stored as (cols_B x rows_B) in row-major order.
 */
export function dotTrans_f64(
    a_ptr: usize,
    b_ptr: usize,
    out_ptr: usize,
    rows_A: i32,
    cols_A: i32,
    cols_B: i32
): void {
    for (let i = 0; i < rows_A; ++i) {
        const a_row = a_ptr + ((i * cols_A) << 3);
        const i_cols_B = i * cols_B;
        for (let j = 0; j < cols_B; ++j) {
            const sum = dot_simd_f64(a_row, b_ptr + ((j * cols_A) << 3), cols_A);
            store<f64>(out_ptr + ((i_cols_B + j) << 3), sum);
        }
    }
}

/**
 * Outer product C (len x len) = A (len) x B (len)
 */
export function outer_f64(a_ptr: usize, b_ptr: usize, out_ptr: usize, len: i32): void {
    for (let i = 0; i < len; ++i) {
        const ai = load<f64>(a_ptr + (i << 3));
        const i_len = i * len;
        for (let j = 0; j < len; ++j) {
            if (i <= j) {
                const bj = load<f64>(b_ptr + (j << 3));
                store<f64>(out_ptr + ((i_len + j) << 3), ai * bj);
            } else {
                const entry_ji = load<f64>(out_ptr + ((j * len + i) << 3));
                store<f64>(out_ptr + ((i_len + j) << 3), entry_ji);
            }
        }
    }
}
