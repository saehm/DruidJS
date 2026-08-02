// src/dimred/SQDMDS.as.ts

import { sqdist_f64, zero_f64 } from "../wasm/shared.as";

/**
 * Computes quartet gradients and updates gradient buffer in WASM memory without GC object allocations.
 */
export function sqdmds_fill_grads_f64(
    y_ptr: usize,           // N x d (f64)
    x_ptr: usize,           // N x D_hd (f64)
    quartets_ptr: usize,    // num_quartets x 4 (i32)
    num_quartets: i32,
    grads_ptr: usize,       // N x d (f64)
    n: i32,
    d_ld: i32,
    d_hd: i32,
    use_exaggeration: boolean,
    is_precomputed: boolean
): void {
    zero_f64(grads_ptr, n * d_ld);
    sqdmds_fill_grads_range_f64(y_ptr, x_ptr, quartets_ptr, num_quartets, grads_ptr, n, d_ld, d_hd, use_exaggeration, is_precomputed, 0, num_quartets);
}

/**
 * High-dimensional distance between rows `a` and `b`.
 */
function sqdmds_hd_distance(
    x_ptr: usize,
    a: i32,
    b: i32,
    n: i32,
    d_hd: i32,
    use_exaggeration: boolean,
    is_precomputed: boolean
): f64 {
    if (is_precomputed) {
        const value = load<f64>(x_ptr + ((a * n + b) << 3));
        return use_exaggeration ? value * value : value;
    }
    const sum = sqdist_f64(x_ptr + ((a * d_hd) << 3), x_ptr + ((b * d_hd) << 3), d_hd);
    // Exaggeration squares the metric, so the square root cancels out.
    return use_exaggeration ? sum : Math.sqrt(sum);
}

/**
 * Gradients for one element of the loss function's sum, accumulated into the four target rows.
 */
function sqdmds_abcd_grads(
    y_ptr: usize,
    grads_ptr: usize,
    a: i32, b: i32, c: i32, d: i32,
    d_ab: f64, d_ac: f64, d_ad: f64, d_bc: f64, d_bd: f64, d_cd: f64,
    p_ab: f64,
    sum_ld: f64,
    dim: i32
): void {
    const ratio = d_ab / sum_ld;
    const twice_ratio = 2.0 * ((p_ab - ratio) / sum_ld);
    const rt = ratio * twice_ratio;

    const a_off = a * dim;
    const b_off = b * dim;
    const c_off = c * dim;
    const d_off = d * dim;

    for (let t = 0; t < dim; ++t) {
        const ya = load<f64>(y_ptr + ((a_off + t) << 3));
        const yb = load<f64>(y_ptr + ((b_off + t) << 3));
        const yc = load<f64>(y_ptr + ((c_off + t) << 3));
        const yd = load<f64>(y_ptr + ((d_off + t) << 3));

        const ab = (ya - yb) / d_ab;
        const ac = (ya - yc) / d_ac;
        const ad = (ya - yd) / d_ad;
        const bc = (yb - yc) / d_bc;
        const bd = (yb - yd) / d_bd;
        const cd = (yc - yd) / d_cd;

        const g_a = ((ab + ac + ad) * ratio - ab) * twice_ratio;
        const g_b = ((-ab + bc + bd) * ratio + ab) * twice_ratio;
        const g_c = (-ac - bc + cd) * rt;
        const g_d = (-ad - bd - cd) * rt;

        store<f64>(grads_ptr + ((a_off + t) << 3), load<f64>(grads_ptr + ((a_off + t) << 3)) + g_a);
        store<f64>(grads_ptr + ((b_off + t) << 3), load<f64>(grads_ptr + ((b_off + t) << 3)) + g_b);
        store<f64>(grads_ptr + ((c_off + t) << 3), load<f64>(grads_ptr + ((c_off + t) << 3)) + g_c);
        store<f64>(grads_ptr + ((d_off + t) << 3), load<f64>(grads_ptr + ((d_off + t) << 3)) + g_d);
    }
}

/**
 * Accumulates quartet gradients over `[start_q, end_q)` without zeroing the buffer first.
 */
export function sqdmds_fill_grads_range_f64(
    y_ptr: usize,
    x_ptr: usize,
    quartets_ptr: usize,
    num_quartets: i32,
    grads_ptr: usize,
    n: i32,
    d_ld: i32,
    d_hd: i32,
    use_exaggeration: boolean,
    is_precomputed: boolean,
    start_q: i32,
    end_q: i32
): void {
    const p = new Array<f64>(6);

    for (let q = start_q; q < end_q; ++q) {
        const q_4 = (q * 4) << 2;
        const i = load<i32>(quartets_ptr + q_4);
        const j = load<i32>(quartets_ptr + q_4 + 4);
        const k = load<i32>(quartets_ptr + q_4 + 8);
        const l = load<i32>(quartets_ptr + q_4 + 12);

        p[0] = sqdmds_hd_distance(x_ptr, i, j, n, d_hd, use_exaggeration, is_precomputed);
        p[1] = sqdmds_hd_distance(x_ptr, i, k, n, d_hd, use_exaggeration, is_precomputed);
        p[2] = sqdmds_hd_distance(x_ptr, i, l, n, d_hd, use_exaggeration, is_precomputed);
        p[3] = sqdmds_hd_distance(x_ptr, j, k, n, d_hd, use_exaggeration, is_precomputed);
        p[4] = sqdmds_hd_distance(x_ptr, j, l, n, d_hd, use_exaggeration, is_precomputed);
        p[5] = sqdmds_hd_distance(x_ptr, k, l, n, d_hd, use_exaggeration, is_precomputed);

        let hd_sum: f64 = 0.0;
        for (let t = 0; t < 6; ++t) hd_sum += p[t];
        if (hd_sum > 0.0) {
            for (let t = 0; t < 6; ++t) {
                p[t] = p[t] / hd_sum + 1e-11;
            }
        }

        // Low-dimensional distances, nudged off zero so the divisions below stay finite.
        const row_i = y_ptr + ((i * d_ld) << 3);
        const row_j = y_ptr + ((j * d_ld) << 3);
        const row_k = y_ptr + ((k * d_ld) << 3);
        const row_l = y_ptr + ((l * d_ld) << 3);
        const d_ij = Math.sqrt(sqdist_f64(row_i, row_j, d_ld)) + 1e-12;
        const d_ik = Math.sqrt(sqdist_f64(row_i, row_k, d_ld)) + 1e-12;
        const d_il = Math.sqrt(sqdist_f64(row_i, row_l, d_ld)) + 1e-12;
        const d_jk = Math.sqrt(sqdist_f64(row_j, row_k, d_ld)) + 1e-12;
        const d_jl = Math.sqrt(sqdist_f64(row_j, row_l, d_ld)) + 1e-12;
        const d_kl = Math.sqrt(sqdist_f64(row_k, row_l, d_ld)) + 1e-12;
        const sum_ld = d_ij + d_ik + d_il + d_jk + d_jl + d_kl;

        // The same gradient for each of the six pairs, with the quartet permuted into place.
        sqdmds_abcd_grads(y_ptr, grads_ptr, i, j, k, l, d_ij, d_ik, d_il, d_jk, d_jl, d_kl, p[0], sum_ld, d_ld);
        sqdmds_abcd_grads(y_ptr, grads_ptr, i, k, j, l, d_ik, d_ij, d_il, d_jk, d_kl, d_jl, p[1], sum_ld, d_ld);
        sqdmds_abcd_grads(y_ptr, grads_ptr, i, l, k, j, d_il, d_ik, d_ij, d_kl, d_jl, d_jk, p[2], sum_ld, d_ld);
        sqdmds_abcd_grads(y_ptr, grads_ptr, j, k, i, l, d_jk, d_ij, d_jl, d_ik, d_kl, d_il, p[3], sum_ld, d_ld);
        sqdmds_abcd_grads(y_ptr, grads_ptr, j, l, i, k, d_jl, d_ij, d_jk, d_il, d_kl, d_ik, p[4], sum_ld, d_ld);
        sqdmds_abcd_grads(y_ptr, grads_ptr, k, l, i, j, d_kl, d_ik, d_jk, d_il, d_jl, d_ij, p[5], sum_ld, d_ld);
    }
}

/**
 * Applies Nesterov momentum update and updates embedding Y in WASM SIMD.
 */
export function sqdmds_nestrov_step_f64(
    y_ptr: usize,          // N x d (f64)
    m_ptr: usize,          // N x d (f64)
    grads_ptr: usize,      // N x d (f64)
    n: i32,
    d: i32,
    lr: f64
): void {
    for (let i = 0; i < n; ++i) {
        const i_d = i * d;

        // The caller has already decayed the momentums by 0.99 in order to evaluate the gradients at
        // the look-ahead position, so this kernel must not decay them a second time.
        let sum_sq: f64 = 0.0;
        for (let k = 0; k < d; ++k) {
            const g_val = load<f64>(grads_ptr + ((i_d + k) << 3));
            sum_sq += g_val * g_val;
        }

        const g_norm = Math.sqrt(sum_sq);
        if (g_norm > 0.0) {
            const mul = lr / g_norm;
            for (let k = 0; k < d; ++k) {
                const idx = (i_d + k) << 3;
                const g_val = load<f64>(grads_ptr + idx);
                store<f64>(m_ptr + idx, load<f64>(m_ptr + idx) - mul * g_val);
            }
        }

        // Y moves by its momentum whether or not this row had a gradient — a row with a zero
        // gradient still carries the momentum accumulated in earlier iterations.
        for (let k = 0; k < d; ++k) {
            const idx = (i_d + k) << 3;
            store<f64>(y_ptr + idx, load<f64>(y_ptr + idx) + load<f64>(m_ptr + idx));
        }
    }
}
