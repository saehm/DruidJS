// src/metrics/manhattan.as.ts

/**
 * Vector Manhattan (L1) Distance SIMD kernel.
 */
export function manhattan_distance_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    let sum: f64 = 0.0;
    let i = 0;
    const len_simd = len - 1;
    let sum_v = f64x2.splat(0.0);

    for (; i < len_simd; i += 2) {
        const a_v = v128.load(a_ptr + (i << 3));
        const b_v = v128.load(b_ptr + (i << 3));
        const diff = f64x2.sub(a_v, b_v);
        const abs_diff = f64x2.abs(diff);
        sum_v = f64x2.add(sum_v, abs_diff);
    }

    sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

    for (; i < len; ++i) {
        const a = load<f64>(a_ptr + (i << 3));
        const b = load<f64>(b_ptr + (i << 3));
        sum += Math.abs(a - b);
    }

    return sum;
}
