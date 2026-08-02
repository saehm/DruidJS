// src/clustering/KMeans.as.ts

import { sqdist_simd_f64 } from "../wasm/shared.as";

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
    for (let i = 0; i < n; ++i) {
        const row_i = x_ptr + ((i * d) << 3);
        let min_dist: f64 = F64.MAX_VALUE;
        let min_k: i32 = 0;

        for (let c = 0; c < k; ++c) {
            // Squared distance throughout: the ordering is the same and the square root is not.
            const sum = sqdist_simd_f64(row_i, c_ptr + ((c * d) << 3), d);
            if (sum < min_dist) {
                min_dist = sum;
                min_k = c;
            }
        }

        store<i32>(assignments_ptr + (i << 2), min_k);
    }
}
