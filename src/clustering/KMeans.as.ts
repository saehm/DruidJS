// src/clustering/KMeans.as.ts

/**
 * k-Means Centroid Assignment SIMD kernel.
 * Assigns each of N data points in X (N x D) to nearest centroid in C (K x D).
 */
export function kmeans_assign_f64(
    x_ptr: usize,
    c_ptr: usize,
    assignments_ptr: usize,
    n: i32,
    k: i32,
    d: i32
): void {
    const d_simd = d - 1;

    for (let i = 0; i < n; ++i) {
        const i_d = i * d;
        let min_dist: f64 = F64.MAX_VALUE;
        let min_k: i32 = 0;

        for (let c = 0; c < k; ++c) {
            const c_d = c * d;
            let sum: f64 = 0.0;
            let j = 0;
            let sum_v = f64x2.splat(0.0);

            for (; j < d_simd; j += 2) {
                const xi_v = v128.load(x_ptr + ((i_d + j) << 3));
                const cj_v = v128.load(c_ptr + ((c_d + j) << 3));
                const diff = f64x2.sub(xi_v, cj_v);
                sum_v = f64x2.add(sum_v, f64x2.mul(diff, diff));
            }

            sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

            for (; j < d; ++j) {
                const diff = load<f64>(x_ptr + ((i_d + j) << 3)) - load<f64>(c_ptr + ((c_d + j) << 3));
                sum += diff * diff;
            }

            if (sum < min_dist) {
                min_dist = sum;
                min_k = c;
            }
        }

        store<i32>(assignments_ptr + (i << 2), min_k);
    }
}
