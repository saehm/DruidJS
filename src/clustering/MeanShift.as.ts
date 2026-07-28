// src/clustering/MeanShift.as.ts

/**
 * WASM SIMD MeanShift Iteration Step Kernel.
 * Computes Gaussian or Flat kernel density shifts for points [start_i, end_i) in linear WASM memory.
 */
export function meanshift_step_range_f64(
    points_ptr: usize,     // N x D (f64)
    out_points_ptr: usize, // N x D (f64)
    n: i32,
    d: i32,
    bandwidth: f64,
    use_gaussian: boolean,
    start_i: i32,
    end_i: i32
): f64 {
    const inv_two_bw_sq = 1.0 / (2.0 * bandwidth * bandwidth);
    const d_simd = d - 1;
    let max_shift: f64 = 0.0;

    for (let i = start_i; i < end_i; ++i) {
        const i_d = i * d;
        let sum_weights: f64 = 0.0;

        // Temporary weighted sum buffer
        const weighted_sum = new Array<f64>(d);
        for (let k = 0; k < d; ++k) weighted_sum[k] = 0.0;

        for (let j = 0; j < n; ++j) {
            const j_d = j * d;
            let sum_sq: f64 = 0.0;
            let k = 0;
            let sum_v = f64x2.splat(0.0);

            for (; k < d_simd; k += 2) {
                const pi_v = v128.load(points_ptr + ((i_d + k) << 3));
                const pj_v = v128.load(points_ptr + ((j_d + k) << 3));
                const diff = f64x2.sub(pi_v, pj_v);
                sum_v = f64x2.add(sum_v, f64x2.mul(diff, diff));
            }
            sum_sq += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

            for (; k < d; ++k) {
                const diff = load<f64>(points_ptr + ((i_d + k) << 3)) - load<f64>(points_ptr + ((j_d + k) << 3));
                sum_sq += diff * diff;
            }

            const dist = Math.sqrt(sum_sq);
            let weight: f64 = 0.0;
            if (use_gaussian) {
                weight = Math.exp(-sum_sq * inv_two_bw_sq);
            } else {
                weight = dist <= bandwidth ? 1.0 : 0.0;
            }

            sum_weights += weight;
            for (let k = 0; k < d; ++k) {
                weighted_sum[k] += weight * load<f64>(points_ptr + ((j_d + k) << 3));
            }
        }

        let shift_sq_sum: f64 = 0.0;
        if (sum_weights == 0.0) {
            for (let k = 0; k < d; ++k) {
                const val = weighted_sum[k];
                shift_sq_sum += val * val;
                store<f64>(out_points_ptr + ((i_d + k) << 3), load<f64>(points_ptr + ((i_d + k) << 3)));
            }
        } else {
            const inv_sum = 1.0 / sum_weights;
            for (let k = 0; k < d; ++k) {
                const curr_x = load<f64>(points_ptr + ((i_d + k) << 3));
                const shift_k = weighted_sum[k] * inv_sum - curr_x;
                shift_sq_sum += shift_k * shift_k;
                store<f64>(out_points_ptr + ((i_d + k) << 3), curr_x + shift_k);
            }
        }

        const shift_norm = Math.sqrt(shift_sq_sum);
        if (shift_norm > max_shift) max_shift = shift_norm;
    }

    return max_shift;
}

export function meanshift_step_f64(
    points_ptr: usize,
    out_points_ptr: usize,
    n: i32,
    d: i32,
    bandwidth: f64,
    use_gaussian: boolean
): f64 {
    return meanshift_step_range_f64(points_ptr, out_points_ptr, n, d, bandwidth, use_gaussian, 0, n);
}
