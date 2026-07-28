// src/numerical/neumair_sum.as.ts

/**
 * SIMD Neumaier Compensated Summation.
 */
export function neumair_sum_f64(v_ptr: usize, len: i32): f64 {
    let sum: f64 = 0.0;
    let c: f64 = 0.0;

    for (let i = 0; i < len; ++i) {
        const val = load<f64>(v_ptr + (i << 3));
        const t = sum + val;
        if (Math.abs(sum) >= Math.abs(val)) {
            c += (sum - t) + val;
        } else {
            c += (val - t) + sum;
        }
        sum = t;
    }

    return sum + c;
}
