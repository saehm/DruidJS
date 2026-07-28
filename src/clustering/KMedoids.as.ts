// src/clustering/KMedoids.as.ts

/**
 * k-Medoids Cluster Assignment Kernel.
 * Assigns each of N points to nearest medoid index based on Distance Matrix D (N x N).
 */
export function kmedoids_assign_f64(
    d_ptr: usize,
    medoids_ptr: usize,
    assignments_ptr: usize,
    n: i32,
    k: i32
): void {
    for (let i = 0; i < n; ++i) {
        const i_n = i * n;
        let min_dist: f64 = F64.MAX_VALUE;
        let min_k: i32 = 0;

        for (let m = 0; m < k; ++m) {
            const medoid_idx = load<i32>(medoids_ptr + (m << 2));
            const dist = load<f64>(d_ptr + ((i_n + medoid_idx) << 3));
            if (dist < min_dist) {
                min_dist = dist;
                min_k = m;
            }
        }

        store<i32>(assignments_ptr + (i << 2), min_k);
    }
}
