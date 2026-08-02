// src/matrix/norm.as.ts

import { sqnorm_simd_f64 } from "../wasm/shared.as";

/**
 * L2 Norm SIMD Kernel.
 */
export function norm_f64(v_ptr: usize, len: i32): f64 {
    return Math.sqrt(sqnorm_simd_f64(v_ptr, len));
}
