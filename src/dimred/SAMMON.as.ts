// src/dimred/SAMMON.as.ts

/**
 * Sammon Mapping Single Iteration Pseudo-Newton Step SIMD Kernel.
 * Computes 1st (e1) and 2nd (e2) stress derivatives and updates embedding coordinates Y (N x D).
 */
export function sammon_step_f64(
    y_ptr: usize,
    d_ptr: usize,
    g_ptr: usize,
    n: i32,
    dim: i32,
    magic: f64
): void {
    const total = n * dim;
    const inv_n = 1.0 / f64(n);

    // Array to store column sums for mean centering
    const sums = new Array<f64>(dim);
    for (let k = 0; k < dim; ++k) sums[k] = 0.0;

    for (let i = 0; i < n; ++i) {
        const i_d = i * dim;
        const i_n = i * n;

        // Temporary buffers for e1 and e2
        const e1 = new Array<f64>(dim);
        const e2 = new Array<f64>(dim);
        for (let k = 0; k < dim; ++k) {
            e1[k] = 0.0;
            e2[k] = 0.0;
        }

        for (let j = 0; j < n; ++j) {
            if (i == j) continue;
            const dX = load<f64>(d_ptr + ((i_n + j) << 3));
            if (dX == 0.0) continue;

            const j_d = j * dim;
            let sum_sq: f64 = 0.0;
            for (let k = 0; k < dim; ++k) {
                const diff = load<f64>(y_ptr + ((i_d + k) << 3)) - load<f64>(y_ptr + ((j_d + k) << 3));
                sum_sq += diff * diff;
            }

            let dY = Math.sqrt(sum_sq);
            if (dY < 1e-6) dY = 1e-6;

            const dq = dX - dY;
            const dr = dX * dY;
            const inv_dr = 1.0 / dr;
            const inv_dY = 1.0 / dY;

            for (let k = 0; k < dim; ++k) {
                const delta_k = load<f64>(y_ptr + ((i_d + k) << 3)) - load<f64>(y_ptr + ((j_d + k) << 3));
                e1[k] += (delta_k * dq) * inv_dr;
                const delta_sq = delta_k * delta_k;
                const term2 = (delta_sq * (1.0 + dq * inv_dY)) * inv_dY;
                e2[k] += (dq - term2) * inv_dr;
            }
        }

        for (let k = 0; k < dim; ++k) {
            const abs_e2 = Math.abs(e2[k]);
            const step = abs_e2 > 0.0 ? (magic * e1[k]) / abs_e2 : 0.0;
            const y_ik = load<f64>(y_ptr + ((i_d + k) << 3));
            const val = y_ik + step;
            store<f64>(g_ptr + ((i_d + k) << 3), val);
            sums[k] += val;
        }
    }

    for (let k = 0; k < dim; ++k) {
        sums[k] *= inv_n;
    }

    // Centering: Y = G - mean
    for (let i = 0; i < n; ++i) {
        const i_d = i * dim;
        for (let k = 0; k < dim; ++k) {
            const val = load<f64>(g_ptr + ((i_d + k) << 3)) - sums[k];
            store<f64>(y_ptr + ((i_d + k) << 3), val);
        }
    }
}
