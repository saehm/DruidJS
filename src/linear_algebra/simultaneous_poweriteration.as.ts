// src/linear_algebra/simultaneous_poweriteration.as.ts

import { dot_simd_f64 } from "../wasm/shared.as";

/**
 * SIMD Matrix-Vector Product Y = A * X (A is rows x cols, X is cols).
 */
export function mat_vec_mul_f64(
    a_ptr: usize,
    x_ptr: usize,
    out_ptr: usize,
    rows: i32,
    cols: i32
): void {
    for (let i = 0; i < rows; ++i) {
        const sum = dot_simd_f64(a_ptr + ((i * cols) << 3), x_ptr, cols);
        store<f64>(out_ptr + (i << 3), sum);
    }
}
