// src/matrix/norm.as.ts

/**
 * L2 Norm SIMD Kernel.
 */
export function norm_f64(v_ptr: usize, len: i32): f64 {
    let sum: f64 = 0.0;
    let i = 0;
    const len_simd = len - 1;
    let sum_v = f64x2.splat(0.0);

    for (; i < len_simd; i += 2) {
        const v = v128.load(v_ptr + (i << 3));
        sum_v = f64x2.add(sum_v, f64x2.mul(v, v));
    }

    sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

    for (; i < len; ++i) {
        const val = load<f64>(v_ptr + (i << 3));
        sum += val * val;
    }

    return Math.sqrt(sum);
}
