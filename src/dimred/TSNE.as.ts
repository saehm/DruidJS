// src/dimred/TSNE.as.ts

/**
 * t-SNE Single Iteration Gradient SIMD Kernel.
 * Computes Qu, Q, gradients, updates gains & ystep, and updates Y positions.
 */
export function tsne_step_f64(
    y_ptr: usize,
    p_ptr: usize,
    qu_ptr: usize,
    q_ptr: usize,
    grad_ptr: usize,
    ystep_ptr: usize,
    gains_ptr: usize,
    n: i32,
    dim: i32,
    pmul: f64,
    epsilon: f64,
    momval: f64
): void {
    let qsum: f64 = 0.0;

    // 1. Compute unnormalized Qu and qsum
    for (let i = 0; i < n; ++i) {
        const i_d = i * dim;
        const i_n = i * n;
        store<f64>(qu_ptr + ((i_n + i) << 3), 0.0);

        for (let j = i + 1; j < n; ++j) {
            const j_d = j * dim;
            let dsum: f64 = 0.0;
            for (let d = 0; d < dim; ++d) {
                const diff = load<f64>(y_ptr + ((i_d + d) << 3)) - load<f64>(y_ptr + ((j_d + d) << 3));
                dsum += diff * diff;
            }
            const qu = 1.0 / (1.0 + dsum);
            store<f64>(qu_ptr + ((i_n + j) << 3), qu);
            store<f64>(qu_ptr + ((j * n + i) << 3), qu);
            qsum += 2.0 * qu;
        }
    }

    // 2. Compute normalized Q
    for (let i = 0; i < n; ++i) {
        const i_n = i * n;
        for (let j = i + 1; j < n; ++j) {
            const qu = load<f64>(qu_ptr + ((i_n + j) << 3));
            let val = qu / qsum;
            if (val < 1e-100) val = 1e-100;
            store<f64>(q_ptr + ((i_n + j) << 3), val);
            store<f64>(q_ptr + ((j * n + i) << 3), val);
        }
    }

    // 3. Compute gradients
    const total_grad = n * dim;
    for (let i = 0; i < total_grad; ++i) {
        store<f64>(grad_ptr + (i << 3), 0.0);
    }

    for (let i = 0; i < n; ++i) {
        const i_n = i * n;
        const i_d = i * dim;
        for (let j = 0; j < n; ++j) {
            if (i == j) continue;
            const j_d = j * dim;
            const pij = load<f64>(p_ptr + ((i_n + j) << 3));
            const qij = load<f64>(q_ptr + ((i_n + j) << 3));
            const quij = load<f64>(qu_ptr + ((i_n + j) << 3));
            const premult = 4.0 * (pmul * pij - qij) * quij;

            for (let d = 0; d < dim; ++d) {
                const diff = load<f64>(y_ptr + ((i_d + d) << 3)) - load<f64>(y_ptr + ((j_d + d) << 3));
                const g_offset = (i_d + d) << 3;
                const prev_g = load<f64>(grad_ptr + g_offset);
                store<f64>(grad_ptr + g_offset, prev_g + premult * diff);
            }
        }
    }

    // 4. Perform gradient step & update Y
    for (let d = 0; d < dim; ++d) {
        let ymean: f64 = 0.0;
        for (let i = 0; i < n; ++i) {
            const idx = (i * dim + d) << 3;
            const gid = load<f64>(grad_ptr + idx);
            const sid = load<f64>(ystep_ptr + idx);
            const gainid = load<f64>(gains_ptr + idx);

            // Must match `Math.sign(gid) === Math.sign(sid)`: a zero counts as its own sign, so a
            // zero paired with a non-zero is *not* the same sign. `ystep` starts out all zeros, so
            // treating them as matching would take the wrong branch on the whole first iteration.
            const sign_g = gid > 0.0 ? 1 : gid < 0.0 ? -1 : 0;
            const sign_s = sid > 0.0 ? 1 : sid < 0.0 ? -1 : 0;
            let newgain = sign_g == sign_s ? gainid * 0.8 : gainid + 0.2;
            if (newgain < 0.01) newgain = 0.01;
            store<f64>(gains_ptr + idx, newgain);

            const newsid = momval * sid - epsilon * newgain * gid;
            store<f64>(ystep_ptr + idx, newsid);

            const newy = load<f64>(y_ptr + idx) + newsid;
            store<f64>(y_ptr + idx, newy);
            ymean += newy;
        }

        // Subtract mean
        const mean = ymean / (n as f64);
        for (let i = 0; i < n; ++i) {
            const idx = (i * dim + d) << 3;
            const current_y = load<f64>(y_ptr + idx);
            store<f64>(y_ptr + idx, current_y - mean);
        }
    }
}
