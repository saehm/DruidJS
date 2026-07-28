// src/metrics/canberra.as.ts

/**
 * Vector Canberra Distance SIMD kernel.
 */
export function canberra_distance_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    let sum: f64 = 0.0;
    for (let i = 0; i < len; ++i) {
        const a = load<f64>(a_ptr + (i << 3));
        const b = load<f64>(b_ptr + (i << 3));
        const num = Math.abs(a - b);
        const denom = Math.abs(a) + Math.abs(b);
        if (denom != 0.0) {
            sum += num / denom;
        }
    }
    return sum;
}
