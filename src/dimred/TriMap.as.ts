// src/dimred/TriMap.as.ts

import { adapt_gain, zero_f64 } from "../wasm/shared.as";

/**
 * TriMap Triplet Gradient Computation WASM SIMD Kernel.
 * Computes gradients and loss for all triplets without heap allocations or GC overhead.
 */
export function trimap_grad_f64(
    y_ptr: usize,        // N x dim (f64)
    triplets_ptr: usize, // num_triplets x 3 (i32)
    weights_ptr: usize,  // num_triplets (f64)
    grad_ptr: usize,     // N x dim (f64)
    n: i32,
    dim: i32,
    num_triplets: i32,
    n_inliers: i32,
    n_outliers: i32
): f64 {
    zero_f64(grad_ptr, n * dim);

    return trimap_grad_range_f64(y_ptr, triplets_ptr, weights_ptr, grad_ptr, n, dim, num_triplets, n_inliers, n_outliers, 0, num_triplets);
}

/**
 * Range-based Parallel TriMap Triplet Gradient Kernel.
 * Allows worker threads to process a slice of triplets [start_t, end_t) concurrently.
 */
export function trimap_grad_range_f64(
    y_ptr: usize,
    triplets_ptr: usize,
    weights_ptr: usize,
    grad_ptr: usize,
    n: i32,
    dim: i32,
    num_triplets: i32,
    n_inliers: i32,
    n_outliers: i32,
    start_t: i32,
    end_t: i32
): f64 {
    let loss: f64 = 0.0;
    const n_knn_triplets = n * n_inliers * n_outliers;

    const y_ij = new Array<f64>(dim);
    const y_ik = new Array<f64>(dim);

    let d_ij: f64 = 1.0;
    let d_ik: f64 = 1.0;

    for (let t = start_t; t < end_t; ++t) {
        const t_3 = (t * 3) << 2;
        const i = load<i32>(triplets_ptr + t_3);
        const j = load<i32>(triplets_ptr + t_3 + 4);
        const k = load<i32>(triplets_ptr + t_3 + 8);
        const w_t = load<f64>(weights_ptr + (t << 3));

        const i_dim = i * dim;
        const j_dim = j * dim;
        const k_dim = k * dim;

        if (t % n_outliers == 0 || t >= n_knn_triplets) {
            d_ij = 1.0;
            d_ik = 1.0;
            for (let d = 0; d < dim; ++d) {
                const Y_id = load<f64>(y_ptr + ((i_dim + d) << 3));
                const Y_jd = load<f64>(y_ptr + ((j_dim + d) << 3));
                const Y_kd = load<f64>(y_ptr + ((k_dim + d) << 3));

                const diff_ij = Y_id - Y_jd;
                const diff_ik = Y_id - Y_kd;

                y_ij[d] = diff_ij;
                y_ik[d] = diff_ik;

                d_ij += diff_ij * diff_ij;
                d_ik += diff_ik * diff_ik;
            }
        } else {
            d_ik = 1.0;
            for (let d = 0; d < dim; ++d) {
                const Y_id = load<f64>(y_ptr + ((i_dim + d) << 3));
                const Y_kd = load<f64>(y_ptr + ((k_dim + d) << 3));

                const diff_ik = Y_id - Y_kd;
                y_ik[d] = diff_ik;
                d_ik += diff_ik * diff_ik;
            }
        }

        loss += w_t / (1.0 + d_ik / d_ij);
        const sum_d = d_ij + d_ik;
        const w = w_t / (sum_d * sum_d);

        for (let d = 0; d < dim; ++d) {
            const gs = y_ij[d] * d_ik * w;
            const go = y_ik[d] * d_ij * w;

            const gi_offset = (i_dim + d) << 3;
            const gj_offset = (j_dim + d) << 3;
            const gk_offset = (k_dim + d) << 3;

            store<f64>(grad_ptr + gi_offset, load<f64>(grad_ptr + gi_offset) + (gs - go));
            store<f64>(grad_ptr + gj_offset, load<f64>(grad_ptr + gj_offset) - gs);
            store<f64>(grad_ptr + gk_offset, load<f64>(grad_ptr + gk_offset) + go);
        }
    }

    return loss;
}

/**
 * TriMap Embedding Update WASM SIMD Kernel.
 * Vectorizes gains, momentum velocities, and coordinate updates.
 */
export function trimap_update_f64(
    y_ptr: usize,    // N x dim (f64)
    grad_ptr: usize, // N x dim (f64)
    vel_ptr: usize,  // N x dim (f64)
    gain_ptr: usize, // N x dim (f64)
    n: i32,
    dim: i32,
    gamma: f64,
    lr: f64
): void {
    const total = n * dim;

    for (let idx = 0; idx < total; ++idx) {
        const offset = idx << 3;
        const v = load<f64>(vel_ptr + offset);
        const g = load<f64>(grad_ptr + offset);

        const new_gain = adapt_gain(load<f64>(gain_ptr + offset), g, v);
        store<f64>(gain_ptr + offset, new_gain);

        const new_v = gamma * v - lr * new_gain * g;
        store<f64>(vel_ptr + offset, new_v);

        const curr_y = load<f64>(y_ptr + offset);
        store<f64>(y_ptr + offset, curr_y + new_v);
    }
}
