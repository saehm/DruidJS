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
 *
 * This is the body of {@link matmul_f64}, which calls it over all rows -- it is not dead code -- but
 * nothing splits it across the worker pool, and the measurements say not to without changing the
 * pool first. Square by square it is worth having: end to end on 8 workers, 1.2x at n=512, 2.4x at
 * n=800, 4.3x at n=1000, and 3.2x to 4.3x once the pool is warm.
 *
 * The shapes this library actually multiplies are the other kind. `SMACOF` and the power iteration
 * do (n x n) * (n x 2), and there parallel measures 0.08x to 0.11x -- an order of magnitude slower.
 * `run_row_range` copies every input into every worker, so an n x n left operand is copied once per
 * worker to produce an n x 2 result, and the copying dwarfs the arithmetic. Note that an operation
 * count cannot separate the two cases: (4000 x 4000) * (4000 x 2) is 32M multiply-accumulates and
 * loses badly, while (256 x 256) * (256 x 256) is 16M and wins. The guard would have to be work per
 * byte copied, and the honest fix is for workers to share one input mapping rather than each taking
 * a private copy.
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
