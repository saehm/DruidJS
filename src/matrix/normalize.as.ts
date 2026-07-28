// src/matrix/normalize.as.ts

/**
 * In-place SIMD Vector Normalization.
 * Normalizes vector v to unit length.
 */
export function normalize_f64(v_ptr: usize, len: i32): void {
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

    const norm_val = Math.sqrt(sum);
    if (norm_val == 0.0) return;

    const inv_norm = 1.0 / norm_val;
    const inv_norm_v = f64x2.splat(inv_norm);

    i = 0;
    for (; i < len_simd; i += 2) {
        const v = v128.load(v_ptr + (i << 3));
        v128.store(v_ptr + (i << 3), f64x2.mul(v, inv_norm_v));
    }

    for (; i < len; ++i) {
        const val = load<f64>(v_ptr + (i << 3));
        store<f64>(v_ptr + (i << 3), val * inv_norm);
    }
}
