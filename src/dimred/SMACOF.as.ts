// src/dimred/SMACOF.as.ts

/**
 * Single iteration SMACOF Guttman Transform SIMD kernel.
 * Computes Z_new = (1/N) * B(Z) * Z in-place in WASM memory.
 */
export function smacof_step_f64(
    z_ptr: usize,             // N x d (f64)
    target_d_ptr: usize,      // N x N (f64)
    z_new_ptr: usize,         // N x d (f64)
    n: i32,
    d: i32
): f64 {
    const inv_n = 1.0 / f64(n);

    // Reset Z_new to 0
    const total_zd = n * d;
    for (let i = 0; i < total_zd; ++i) {
        store<f64>(z_new_ptr + (i << 3), 0.0);
    }

    for (let i = 0; i < n; ++i) {
        const i_d = i * d;
        const i_n = i * n;
        let bii: f64 = 0.0;

        for (let j = 0; j < n; ++j) {
            if (i == j) continue;
            const j_d = j * d;

            // Compute distance in embedding space dist_Z
            let sum_sq: f64 = 0.0;
            for (let k = 0; k < d; ++k) {
                const diff = load<f64>(z_ptr + ((i_d + k) << 3)) - load<f64>(z_ptr + ((j_d + k) << 3));
                sum_sq += diff * diff;
            }
            const dist_z = Math.sqrt(sum_sq);
            const dist_target = load<f64>(target_d_ptr + ((i_n + j) << 3));

            let bij: f64 = 0.0;
            if (dist_z > 1e-12) {
                bij = -dist_target / dist_z;
            }
            bii -= bij;

            // Multiply row B[i, j] directly into Z_new[i, k]
            for (let k = 0; k < d; ++k) {
                const z_jk = load<f64>(z_ptr + ((j_d + k) << 3));
                const current = load<f64>(z_new_ptr + ((i_d + k) << 3));
                store<f64>(z_new_ptr + ((i_d + k) << 3), current + bij * z_jk);
            }
        }

        // Add diagonal B[i, i] * Z[i, k]
        for (let k = 0; k < d; ++k) {
            const z_ik = load<f64>(z_ptr + ((i_d + k) << 3));
            const current = load<f64>(z_new_ptr + ((i_d + k) << 3));
            const val = (current + bii * z_ik) * inv_n;
            store<f64>(z_new_ptr + ((i_d + k) << 3), val);
        }
    }

    // Stress is reported for the embedding this step produced, not the one it started from — the
    // caller compares it against the previous iteration's stress to decide when to stop, so
    // measuring the pre-step Z here would lag the convergence check by one iteration.
    let stress_num: f64 = 0.0;
    let stress_den: f64 = 0.0;

    for (let i = 0; i < n; ++i) {
        const i_d = i * d;
        const i_n = i * n;
        for (let j = i + 1; j < n; ++j) {
            const j_d = j * d;
            let sum_sq: f64 = 0.0;
            for (let k = 0; k < d; ++k) {
                const diff = load<f64>(z_new_ptr + ((i_d + k) << 3)) - load<f64>(z_new_ptr + ((j_d + k) << 3));
                sum_sq += diff * diff;
            }
            const dist_target = load<f64>(target_d_ptr + ((i_n + j) << 3));
            const diff_target = dist_target - Math.sqrt(sum_sq);
            stress_num += diff_target * diff_target;
            stress_den += dist_target * dist_target;
        }
    }

    if (stress_den < 1e-12) stress_den = 1e-12;
    return Math.sqrt(stress_num / stress_den);
}
