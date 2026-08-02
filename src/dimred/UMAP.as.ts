// src/dimred/UMAP.as.ts

import { sqdist_f64 } from "../wasm/shared.as";

/** Matches `UMAP._clip`. */
function clip(x: f64): f64 {
    if (x > 4.0) return 4.0;
    if (x < -4.0) return -4.0;
    return x;
}

/**
 * Moves the two rows of a pair apart along their connecting axis, `cur` forwards and `oth` back.
 *
 * The attractive and repulsive passes differ only in the sign carried by `grad_coeff`, so both go
 * through here.
 */
function apply_pair_gradient(cur_off: usize, oth_off: usize, grad_coeff: f64, alpha: f64, dim: i32): void {
    for (let d = 0; d < dim; ++d) {
        const c_idx = cur_off + (d << 3);
        const o_idx = oth_off + (d << 3);
        const grad_d = clip(grad_coeff * (load<f64>(c_idx) - load<f64>(o_idx))) * alpha;
        store<f64>(c_idx, load<f64>(c_idx) + grad_d);
        store<f64>(o_idx, load<f64>(o_idx) - grad_d);
    }
}

/**
 * UMAP Single Epoch SGD Optimization Kernel.
 *
 * Mirrors `UMAP._optimize_layout` exactly, including the attractive pass over each active edge and
 * the repulsive pass over its negative samples.
 *
 * Negative samples are *not* drawn here. The JS side owns the seeded random number generator, so it
 * draws them in edge order beforehand and passes the resolved node indices in `neg_ptr`; this kernel
 * walks them with a cursor. Recomputing the per-edge sample count from the same f32 inputs is
 * deterministic, so the cursor stays in step with the pre-pass.
 *
 * The epoch bookkeeping arrays are f32 because the JS side holds them in `Float32Array`s — using f64
 * here would accumulate differently and drift away from the fallback path.
 *
 * `head_embedding` and `tail_embedding` are the same buffer in UMAP, so updates alias deliberately.
 */
export function umap_optimize_epoch_f64(
    embedding_ptr: usize,   // n_points x dim (f64)
    head_ptr: usize,        // n_edges (i32)
    tail_ptr: usize,        // n_edges (i32)
    eps_ptr: usize,         // n_edges (f32) epochs_per_sample
    eons_ptr: usize,        // n_edges (f32) epoch_of_next_sample
    epns_ptr: usize,        // n_edges (f32) epochs_per_negative_sample
    eonns_ptr: usize,       // n_edges (f32) epoch_of_next_negative_sample
    neg_ptr: usize,         // n_neg (i32) pre-drawn negative sample node indices
    n_edges: i32,
    n_neg: i32,
    dim: i32,
    iter: f64,
    a: f64,
    b: f64,
    gamma: f64,
    alpha: f64
): void {
    let cursor = 0;

    for (let i = 0; i < n_edges; ++i) {
        if (f64(load<f32>(eons_ptr + (i << 2))) > iter) continue;

        const j = load<i32>(head_ptr + (i << 2));
        const k = load<i32>(tail_ptr + (i << 2));
        const cur_off = embedding_ptr + ((j * dim) << 3);
        const oth_off = embedding_ptr + ((k * dim) << 3);

        // --- attractive force ---
        const dist = sqdist_f64(cur_off, oth_off, dim);

        if (dist > 0.0) {
            const grad_coeff = (-2.0 * a * b * Math.pow(dist, b - 1.0)) / (a * Math.pow(dist, b) + 1.0);
            apply_pair_gradient(cur_off, oth_off, grad_coeff, alpha, dim);
        }

        const eps = load<f32>(eps_ptr + (i << 2));
        store<f32>(eons_ptr + (i << 2), load<f32>(eons_ptr + (i << 2)) + eps);

        // --- repulsive force over the pre-drawn negative samples ---
        const epns = f64(load<f32>(epns_ptr + (i << 2)));
        const n_neg_samples = (iter - f64(load<f32>(eonns_ptr + (i << 2)))) / epns;

        for (let p: f64 = 0.0; p < n_neg_samples; p += 1.0) {
            if (cursor >= n_neg) break; // pre-pass and kernel disagreed; stop rather than read garbage
            const other = load<i32>(neg_ptr + (cursor << 2));
            ++cursor;
            const neg_off = embedding_ptr + ((other * dim) << 3);

            const ndist = sqdist_f64(cur_off, neg_off, dim);

            if (ndist > 0.0) {
                const grad_coeff = (2.0 * gamma * b) / ((0.01 + ndist) * (a * Math.pow(ndist, b) + 1.0));
                apply_pair_gradient(cur_off, neg_off, grad_coeff, alpha, dim);
            }
        }

        store<f32>(eonns_ptr + (i << 2), f32(f64(load<f32>(eonns_ptr + (i << 2))) + n_neg_samples * epns));
    }
}
