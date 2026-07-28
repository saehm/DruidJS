// src/matrix/distance_matrix.as.ts

/**
 * Pairwise Euclidean distance matrix computation for X (n x d).
 * Computes symmetric n x n distance matrix D.
 *
 * Only the upper triangle is computed; each distance is mirrored into the lower triangle, which
 * halves the work compared to visiting every ordered pair.
 */
export function euclidean_distance_matrix_f64(
    x_ptr: usize,
    out_ptr: usize,
    n: i32,
    d: i32
): void {
    const d_simd = d - 1;

    for (let i = 0; i < n; ++i) {
        const i_d = i * d;
        const i_n = i * n;

        store<f64>(out_ptr + ((i_n + i) << 3), 0.0);

        for (let j = i + 1; j < n; ++j) {
            const j_d = j * d;
            let sum: f64 = 0.0;
            let k = 0;

            let sum_v = f64x2.splat(0.0);

            for (; k < d_simd; k += 2) {
                const xi_v = v128.load(x_ptr + ((i_d + k) << 3));
                const xj_v = v128.load(x_ptr + ((j_d + k) << 3));
                const diff = f64x2.sub(xi_v, xj_v);
                const diff_sq = f64x2.mul(diff, diff);
                sum_v = f64x2.add(sum_v, diff_sq);
            }

            sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

            for (; k < d; ++k) {
                const diff = load<f64>(x_ptr + ((i_d + k) << 3)) - load<f64>(x_ptr + ((j_d + k) << 3));
                sum += diff * diff;
            }

            const dist = Math.sqrt(sum);
            store<f64>(out_ptr + ((i_n + j) << 3), dist);
            store<f64>(out_ptr + ((j * n + i) << 3), dist);
        }
    }
}

/**
 * Multi-threaded Range-based Pairwise Euclidean Distance Matrix Kernel.
 * Allows worker threads to compute distance matrix rows [start_row, end_row) in parallel.
 */
export function euclidean_distance_matrix_range_f64(
    x_ptr: usize,
    out_ptr: usize,
    n: i32,
    d: i32,
    start_row: i32,
    end_row: i32
): void {
    const d_simd = d - 1;

    for (let i = start_row; i < end_row; ++i) {
        const i_d = i * d;
        const i_n = i * n;

        store<f64>(out_ptr + ((i_n + i) << 3), 0.0);

        for (let j = 0; j < n; ++j) {
            if (i == j) continue;
            const j_d = j * d;
            let sum: f64 = 0.0;
            let k = 0;

            let sum_v = f64x2.splat(0.0);

            for (; k < d_simd; k += 2) {
                const xi_v = v128.load(x_ptr + ((i_d + k) << 3));
                const xj_v = v128.load(x_ptr + ((j_d + k) << 3));
                const diff = f64x2.sub(xi_v, xj_v);
                const diff_sq = f64x2.mul(diff, diff);
                sum_v = f64x2.add(sum_v, diff_sq);
            }

            sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

            for (; k < d; ++k) {
                const diff = load<f64>(x_ptr + ((i_d + k) << 3)) - load<f64>(x_ptr + ((j_d + k) << 3));
                sum += diff * diff;
            }

            const dist = Math.sqrt(sum);
            store<f64>(out_ptr + ((i_n + j) << 3), dist);
        }
    }
}
