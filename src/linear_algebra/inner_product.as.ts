// src/linear_algebra/inner_product.as.ts

import { dot_simd_f64 } from "../wasm/shared.as";

/**
 * SIMD Vector Inner Product: <a, b> = sum(a_i * b_i).
 */
export function inner_product_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    return dot_simd_f64(a_ptr, b_ptr, len);
}
