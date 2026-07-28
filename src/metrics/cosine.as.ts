// src/metrics/cosine.as.ts

/**
 * Vector Cosine Distance SIMD kernel.
 */
export function cosine_distance_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    let dot: f64 = 0.0;
    let normA: f64 = 0.0;
    let normB: f64 = 0.0;

    let i = 0;
    const len_simd = len - 1;
    let dot_v = f64x2.splat(0.0);
    let normA_v = f64x2.splat(0.0);
    let normB_v = f64x2.splat(0.0);

    for (; i < len_simd; i += 2) {
        const a_v = v128.load(a_ptr + (i << 3));
        const b_v = v128.load(b_ptr + (i << 3));

        dot_v = f64x2.add(dot_v, f64x2.mul(a_v, b_v));
        normA_v = f64x2.add(normA_v, f64x2.mul(a_v, a_v));
        normB_v = f64x2.add(normB_v, f64x2.mul(b_v, b_v));
    }

    dot += f64x2.extract_lane(dot_v, 0) + f64x2.extract_lane(dot_v, 1);
    normA += f64x2.extract_lane(normA_v, 0) + f64x2.extract_lane(normA_v, 1);
    normB += f64x2.extract_lane(normB_v, 0) + f64x2.extract_lane(normB_v, 1);

    for (; i < len; ++i) {
        const a = load<f64>(a_ptr + (i << 3));
        const b = load<f64>(b_ptr + (i << 3));
        dot += a * b;
        normA += a * a;
        normB += b * b;
    }

    // Signals "no direction" to the JS wrapper, which reports such pairs as orthogonal.
    if (normA == 0.0 || normB == 0.0) return 1.0;
    // One square root of the product rather than the product of two square roots: the latter rounds
    // twice, so identical vectors would not come out at a similarity of exactly 1.
    return 1.0 - (dot / Math.sqrt(normA * normB));
}
