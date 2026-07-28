// src/linear_algebra/inner_product.as.ts

/**
 * SIMD Vector Inner Product: <a, b> = sum(a_i * b_i).
 */
export function inner_product_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    let sum: f64 = 0.0;
    let i = 0;
    const len_simd = len - 1;
    let sum_v = f64x2.splat(0.0);

    for (; i < len_simd; i += 2) {
        const a_v = v128.load(a_ptr + (i << 3));
        const b_v = v128.load(b_ptr + (i << 3));
        sum_v = f64x2.add(sum_v, f64x2.mul(a_v, b_v));
    }

    sum += f64x2.extract_lane(sum_v, 0) + f64x2.extract_lane(sum_v, 1);

    for (; i < len; ++i) {
        const a = load<f64>(a_ptr + (i << 3));
        const b = load<f64>(b_ptr + (i << 3));
        sum += a * b;
    }

    return sum;
}
