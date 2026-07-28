// src/metrics/chebyshev.as.ts

/**
 * Vector Chebyshev (L_infinity) Distance SIMD kernel.
 */
export function chebyshev_distance_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    let max_diff: f64 = 0.0;
    let i = 0;
    const len_simd = len - 1;

    for (; i < len_simd; i += 2) {
        const a_v = v128.load(a_ptr + (i << 3));
        const b_v = v128.load(b_ptr + (i << 3));
        const diff = f64x2.abs(f64x2.sub(a_v, b_v));
        const d0 = f64x2.extract_lane(diff, 0);
        const d1 = f64x2.extract_lane(diff, 1);
        if (d0 > max_diff) max_diff = d0;
        if (d1 > max_diff) max_diff = d1;
    }

    for (; i < len; ++i) {
        const diff = Math.abs(load<f64>(a_ptr + (i << 3)) - load<f64>(b_ptr + (i << 3)));
        if (diff > max_diff) max_diff = diff;
    }

    return max_diff;
}
