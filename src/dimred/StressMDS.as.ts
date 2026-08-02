// src/dimred/StressMDS.as.ts

import { sqdist_f64, zero_f64 } from "../wasm/shared.as";

/**
 * Weighted stress of the configuration at `cfg_ptr`, `sum_{i<j} w_ij (||y_i - y_j|| - d_ij)^2`.
 *
 * The line search in {@link stress_mds_step_f64} accepts a trial step by comparing its energy
 * against the incumbent's, so the two must be measured the same way down to the summation order.
 * Keeping one definition is what makes that structural rather than a thing to remember.
 */
function weighted_energy(
    cfg_ptr: usize,
    target_d_ptr: usize,
    w_ptr: usize,
    n: i32,
    dim: i32
): f64 {
    let energy: f64 = 0.0;
    for (let i = 0; i < n; ++i) {
        const i_dim = i * dim;
        const i_n = i * n;
        for (let j = i + 1; j < n; ++j) {
            const w_ij = load<f64>(w_ptr + ((i_n + j) << 3));
            if (!(w_ij > 0.0)) continue;
            const d_ij = load<f64>(target_d_ptr + ((i_n + j) << 3));
            const dist = Math.sqrt(sqdist_f64(cfg_ptr + (i_dim << 3), cfg_ptr + ((j * dim) << 3), dim));
            const residual = dist - d_ij;
            energy += w_ij * residual * residual;
        }
    }
    return energy;
}

/**
 * One weighted-stress iteration: preconditioned gradient, backtracking line search, and the update.
 *
 * Weights arrive precomputed as an N x N matrix rather than as an exponent evaluated here. That
 * keeps `Math.pow` out of the inner loop, and it is what lets the same kernel serve an exponent, an
 * explicit weight matrix (with zeros marking missing pairs) and a caller-supplied function without
 * any of them being a special case. A pair is skipped whenever its weight is not positive.
 *
 * The whole line search lives in the kernel rather than just the gradient, because the search
 * evaluates the energy of a trial configuration up to twelve times per iteration. Returning the
 * gradient to JS and stepping there would copy the N x dim embedding across the boundary once per
 * trial; keeping the loop here copies it once per accepted step.
 *
 * `y_ptr` is updated in place when a step is accepted. `out_ptr` receives two f64s: the adapted step
 * size, and 1.0 if the step was accepted (0.0 if the search ran out of trials or the gradient
 * vanished, either of which means the caller should stop).
 *
 * Returns the energy of the configuration `y_ptr` holds on exit.
 */
export function stress_mds_step_f64(
    y_ptr: usize,             // N x dim (f64), updated in place
    target_d_ptr: usize,      // N x N (f64)
    w_ptr: usize,             // N x N (f64) weights
    g_ptr: usize,             // N x dim (f64) gradient scratch
    trial_ptr: usize,         // N x dim (f64) trial scratch
    n: i32,
    dim: i32,
    step: f64,
    energy_current: f64,
    out_ptr: usize
): f64 {
    const total = n * dim;

    zero_f64(g_ptr, total);

    // Gradient of sum_{i<j} w_ij (||y_i - y_j|| - d_ij)^2, divided per point by its weighted degree
    // sum_j w_ij (Jacobi preconditioning).
    //
    // The division is a units fix, not a tuning trick: it is what makes `step` dimensionless and
    // comparable across weightings. Without it one global step size has to serve points whose
    // weighted degrees differ by orders of magnitude -- what a strongly negative exponent produces
    // between dense and sparse regions -- and the search stalls. Measured over 90 problems at
    // exponent -2, dropping it raised the mean stress ratio from 1.002 to 1.109.
    for (let i = 0; i < n; ++i) {
        const i_dim = i * dim;
        const i_n = i * n;
        const row_i = y_ptr + (i_dim << 3);
        let degree: f64 = 0.0;
        for (let j = 0; j < n; ++j) {
            if (i == j) continue;
            const w_ij = load<f64>(w_ptr + ((i_n + j) << 3));
            if (!(w_ij > 0.0)) continue;
            degree += w_ij;

            const d_ij = load<f64>(target_d_ptr + ((i_n + j) << 3));
            const row_j = y_ptr + ((j * dim) << 3);
            const dist = Math.sqrt(sqdist_f64(row_i, row_j, dim));
            // Coincident embedded points carry no direction; the other pairs separate them next step.
            if (dist < 1e-12) continue;

            const coeff = 2.0 * w_ij * (dist - d_ij) / dist;
            for (let k = 0; k < dim; ++k) {
                const delta = load<f64>(row_i + (k << 3)) - load<f64>(row_j + (k << 3));
                const current = load<f64>(g_ptr + ((i_dim + k) << 3));
                store<f64>(g_ptr + ((i_dim + k) << 3), current + coeff * delta);
            }
        }

        if (degree > 0.0) {
            const inv_degree = 1.0 / degree;
            for (let k = 0; k < dim; ++k) {
                store<f64>(g_ptr + ((i_dim + k) << 3), load<f64>(g_ptr + ((i_dim + k) << 3)) * inv_degree);
            }
        }
    }

    let g_norm: f64 = 0.0;
    for (let i = 0; i < total; ++i) {
        const g = load<f64>(g_ptr + (i << 3));
        g_norm += g * g;
    }
    g_norm = Math.sqrt(g_norm);

    if (g_norm < 1e-12) {
        store<f64>(out_ptr, step);
        store<f64>(out_ptr + 8, 0.0);
        return energy_current;
    }

    for (let trial = 0; trial < 12; ++trial) {
        // The preconditioned gradient is already in configuration units, so `step` multiplies it
        // directly rather than being normalised by the gradient norm.
        for (let i = 0; i < total; ++i) {
            store<f64>(
                trial_ptr + (i << 3),
                load<f64>(y_ptr + (i << 3)) - step * load<f64>(g_ptr + (i << 3))
            );
        }

        const energy = weighted_energy(trial_ptr, target_d_ptr, w_ptr, n, dim);

        if (energy < energy_current) {
            for (let i = 0; i < total; ++i) {
                store<f64>(y_ptr + (i << 3), load<f64>(trial_ptr + (i << 3)));
            }
            store<f64>(out_ptr, step * 1.2);
            store<f64>(out_ptr + 8, 1.0);
            return energy;
        }
        step *= 0.5;
    }

    store<f64>(out_ptr, step);
    store<f64>(out_ptr + 8, 0.0);
    return energy_current;
}

/**
 * Weighted stress of a configuration, `sum_{i<j} w_ij (||y_i - y_j|| - d_ij)^2`.
 *
 * Needed once before the iteration loop to seed `energy_current`; every later value comes back from
 * `stress_mds_step_f64` directly.
 */
export function stress_mds_energy_f64(
    y_ptr: usize,
    target_d_ptr: usize,
    w_ptr: usize,
    n: i32,
    dim: i32
): f64 {
    return weighted_energy(y_ptr, target_d_ptr, w_ptr, n, dim);
}
