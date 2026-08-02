// src/matrix/normalize.as.ts

import { sqnorm_simd_f64 } from "../wasm/shared.as";

/**
 * In-place SIMD Vector Normalization.
 * Normalizes vector v to unit length.
 */
export function normalize_f64(v_ptr: usize, len: i32): void {
    const norm_val = Math.sqrt(sqnorm_simd_f64(v_ptr, len));
    if (norm_val == 0.0) return;

    const inv_norm = 1.0 / norm_val;
    const inv_norm_v = f64x2.splat(inv_norm);
    const len_simd = len - 1;

    let i = 0;
    for (; i < len_simd; i += 2) {
        const v = v128.load(v_ptr + (i << 3));
        v128.store(v_ptr + (i << 3), f64x2.mul(v, inv_norm_v));
    }

    for (; i < len; ++i) {
        const val = load<f64>(v_ptr + (i << 3));
        store<f64>(v_ptr + (i << 3), val * inv_norm);
    }
}
