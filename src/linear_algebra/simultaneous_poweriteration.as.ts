// src/linear_algebra/simultaneous_poweriteration.as.ts

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
    const cols_simd = cols - 1;

    for (let i = 0; i < rows; ++i) {
        const i_cols = i * cols;
        let sum: f64 = 0.0;
        let j = 0;
        let sum_v = f64x2.splat(0.0);

        for (; j < cols_simd; j += 2) {
            const a_v = v128.load(a_ptr + ((i_cols + j) << 3));
            const x_v = v128.load(x_ptr + (j << 3));
            sum_v = f64x2.add(sum_v, f64x2.mul(a_v, x_v));
        }

        sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

        for (; j < cols; ++j) {
            const a = load<f64>(a_ptr + ((i_cols + j) << 3));
            const x = load<f64>(x_ptr + (j << 3));
            sum += a * x;
        }

        store<f64>(out_ptr + (i << 3), sum);
    }
}
