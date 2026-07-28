// src/metrics/bray_curtis.as.ts

/**
 * Vector Bray-Curtis Distance SIMD kernel.
 */
export function bray_curtis_distance_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    let num_sum: f64 = 0.0;
    let denom_sum: f64 = 0.0;

    let i = 0;
    const len_simd = len - 1;
    let num_v = f64x2.splat(0.0);
    let denom_v = f64x2.splat(0.0);

    for (; i < len_simd; i += 2) {
        const a_v = v128.load(a_ptr + (i << 3));
        const b_v = v128.load(b_ptr + (i << 3));
        const diff_abs = f64x2.abs(f64x2.sub(a_v, b_v));
        const sum_abs = f64x2.abs(f64x2.add(a_v, b_v));

        num_v = f64x2.add(num_v, diff_abs);
        denom_v = f64x2.add(denom_v, sum_abs);
    }

    num_sum += f64x2.extract_lane(num_v, 0) + f64x2.extract_lane(num_v, 1);
    denom_sum += f64x2.extract_lane(denom_v, 0) + f64x2.extract_lane(denom_v, 1);

    for (; i < len; ++i) {
        const a = load<f64>(a_ptr + (i << 3));
        const b = load<f64>(b_ptr + (i << 3));
        num_sum += Math.abs(a - b);
        denom_sum += Math.abs(a + b);
    }

    if (denom_sum == 0.0) return 0.0;
    return num_sum / denom_sum;
}
