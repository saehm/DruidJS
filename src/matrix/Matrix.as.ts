// src/matrix/Matrix.as.ts

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
    const cols_B_simd = cols_B - 1;

    for (let i = start_row; i < end_row; ++i) {
        const i_cols_A = i * cols_A;
        const i_cols_B = i * cols_B;

        for (let j = 0; j < cols_B; ++j) {
            store<f64>(out_ptr + ((i_cols_B + j) << 3), 0.0);
        }

        for (let k = 0; k < cols_A; ++k) {
            const aik = load<f64>(a_ptr + ((i_cols_A + k) << 3));
            if (aik == 0.0) continue;

            const k_cols_B = k * cols_B;
            const aik_v = f64x2.splat(aik);

            let j = 0;
            for (; j < cols_B_simd; j += 2) {
                const b_idx = (k_cols_B + j) << 3;
                const out_idx = (i_cols_B + j) << 3;

                const b_v = v128.load(b_ptr + b_idx);
                const out_v = v128.load(out_ptr + out_idx);

                const prod = f64x2.mul(aik_v, b_v);
                const res = f64x2.add(out_v, prod);

                v128.store(out_ptr + out_idx, res);
            }

            for (; j < cols_B; ++j) {
                const bkj = load<f64>(b_ptr + ((k_cols_B + j) << 3));
                const out_offset = (i_cols_B + j) << 3;
                const prev = load<f64>(out_ptr + out_offset);
                store<f64>(out_ptr + out_offset, prev + aik * bkj);
            }
        }
    }
}

/**
 * Transposed Matrix Multiplication: C (rows_A x cols_B) = A^T * B
 * A is stored as (cols_A x rows_A) in row-major order.
 */
export function transDot_f64(
    a_ptr: usize,
    b_ptr: usize,
    out_ptr: usize,
    cols_A: i32,
    rows_A: i32,
    cols_B: i32
): void {
    const total_out = rows_A * cols_B;
    for (let i = 0; i < total_out; ++i) {
        store<f64>(out_ptr + (i << 3), 0.0);
    }

    for (let k = 0; k < cols_A; ++k) {
        const k_rows_A = k * rows_A;
        const k_cols_B = k * cols_B;
        for (let i = 0; i < rows_A; ++i) {
            const aki = load<f64>(a_ptr + ((k_rows_A + i) << 3));
            if (aki == 0.0) continue;
            const aki_v = f64x2.splat(aki);
            const cols_B_simd = cols_B - 1;

            let j = 0;
            for (; j < cols_B_simd; j += 2) {
                const b_idx = (k_cols_B + j) << 3;
                const out_idx = (i * cols_B + j) << 3;
                const b_v = v128.load(b_ptr + b_idx);
                const out_v = v128.load(out_ptr + out_idx);
                v128.store(out_ptr + out_idx, f64x2.add(out_v, f64x2.mul(aki_v, b_v)));
            }

            for (; j < cols_B; ++j) {
                const bkj = load<f64>(b_ptr + ((k_cols_B + j) << 3));
                const out_offset = (i * cols_B + j) << 3;
                const prev = load<f64>(out_ptr + out_offset);
                store<f64>(out_ptr + out_offset, prev + aki * bkj);
            }
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
    const cols_A_simd = cols_A - 1;

    for (let i = 0; i < rows_A; ++i) {
        const i_cols_A = i * cols_A;
        const i_cols_B = i * cols_B;
        for (let j = 0; j < cols_B; ++j) {
            const j_cols_A = j * cols_A;
            let sum: f64 = 0.0;
            let k = 0;
            let sum_v = f64x2.splat(0.0);

            for (; k < cols_A_simd; k += 2) {
                const a_v = v128.load(a_ptr + ((i_cols_A + k) << 3));
                const b_v = v128.load(b_ptr + ((j_cols_A + k) << 3));
                sum_v = f64x2.add(sum_v, f64x2.mul(a_v, b_v));
            }

            sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

            for (; k < cols_A; ++k) {
                const aik = load<f64>(a_ptr + ((i_cols_A + k) << 3));
                const bjk = load<f64>(b_ptr + ((j_cols_A + k) << 3));
                sum += aik * bjk;
            }

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
