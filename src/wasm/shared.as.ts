// src/wasm/shared.as.ts

/**
 * Building blocks shared by the kernels.
 *
 * Deliberately not listed in `entry.as.ts`: nothing here is a kernel, and re-exporting would put a
 * call boundary in the middle of an inner loop.
 *
 * These carried an explicit `@inline` at first, on the assumption that a call in an inner loop had
 * to be paid for. Measured, that was wrong twice over: `asc -O3` already inlines bodies this small,
 * so forcing it changed no kernel's timing beyond run-to-run noise, and it cost 1037 bytes of
 * duplicated code -- enough to make a refactor that deletes 400 lines of source grow the binary.
 * Output is bit-identical either way. Leave the decision to the optimiser unless a benchmark says
 * otherwise.
 */

/**
 * Horizontal sum of an `f64x2` accumulator, low lane first.
 *
 * Kept `@inline` where the others are not: a `v128` parameter draws AS112 ("exchange of 'v128'
 * values is not supported by all embeddings") as soon as this compiles to a real function. It is
 * harmless here -- the function is not in `entry.as.ts`, so no `v128` ever crosses the host boundary
 * -- but inlining two lane extracts costs nothing and keeps the build warning-free.
 */
@inline
export function hsum_f64x2(v: v128): f64 {
    return f64x2.extract_lane(v, 0) + f64x2.extract_lane(v, 1);
}

/**
 * Squared Euclidean distance between the `len`-long rows at `a_ptr` and `b_ptr`, accumulated in
 * SIMD lanes.
 *
 * Pairs whose lanes are combined only at the end, so for `len >= 4` the result differs in the last
 * bit from {@link sqdist_f64}. Use this wherever the loop is over the *input* dimensionality, which
 * is where the vectorisation pays for itself.
 */
export function sqdist_simd_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    let sum: f64 = 0.0;
    let i = 0;
    const len_simd = len - 1;
    let sum_v = f64x2.splat(0.0);

    for (; i < len_simd; i += 2) {
        const a_v = v128.load(a_ptr + (i << 3));
        const b_v = v128.load(b_ptr + (i << 3));
        const diff = f64x2.sub(a_v, b_v);
        sum_v = f64x2.add(sum_v, f64x2.mul(diff, diff));
    }

    sum += hsum_f64x2(sum_v);

    for (; i < len; ++i) {
        const diff = load<f64>(a_ptr + (i << 3)) - load<f64>(b_ptr + (i << 3));
        sum += diff * diff;
    }

    return sum;
}

/**
 * Squared Euclidean distance between the `len`-long rows at `a_ptr` and `b_ptr`, accumulated
 * left to right.
 *
 * The optimisers use this rather than {@link sqdist_simd_f64} because their JS fallbacks sum the
 * same way, and the two paths are compared against each other. Up to `len == 3` the orders coincide
 * and the choice is free; above that they do not, and `UMAP` and `SAMMON` are chaotic enough that a
 * last-bit difference grows into a visibly different layout — see the reproducibility note in
 * `src/wasm/index.js`.
 */
export function sqdist_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    let sum: f64 = 0.0;
    for (let i = 0; i < len; ++i) {
        const diff = load<f64>(a_ptr + (i << 3)) - load<f64>(b_ptr + (i << 3));
        sum += diff * diff;
    }
    return sum;
}

/** Sum of squares of the `len`-long vector at `v_ptr`, accumulated in SIMD lanes. */
export function sqnorm_simd_f64(v_ptr: usize, len: i32): f64 {
    let sum: f64 = 0.0;
    let i = 0;
    const len_simd = len - 1;
    let sum_v = f64x2.splat(0.0);

    for (; i < len_simd; i += 2) {
        const v = v128.load(v_ptr + (i << 3));
        sum_v = f64x2.add(sum_v, f64x2.mul(v, v));
    }

    sum += hsum_f64x2(sum_v);

    for (; i < len; ++i) {
        const val = load<f64>(v_ptr + (i << 3));
        sum += val * val;
    }

    return sum;
}

/** Inner product of the `len`-long vectors at `a_ptr` and `b_ptr`, accumulated in SIMD lanes. */
export function dot_simd_f64(a_ptr: usize, b_ptr: usize, len: i32): f64 {
    let sum: f64 = 0.0;
    let i = 0;
    const len_simd = len - 1;
    let sum_v = f64x2.splat(0.0);

    for (; i < len_simd; i += 2) {
        const a_v = v128.load(a_ptr + (i << 3));
        const b_v = v128.load(b_ptr + (i << 3));
        sum_v = f64x2.add(sum_v, f64x2.mul(a_v, b_v));
    }

    sum += hsum_f64x2(sum_v);

    for (; i < len; ++i) {
        const a = load<f64>(a_ptr + (i << 3));
        const b = load<f64>(b_ptr + (i << 3));
        sum += a * b;
    }

    return sum;
}

/**
 * `out[0..len) += scale * b[0..len)`, accumulated in SIMD lanes.
 *
 * The inner two loops of every "accumulate a scaled row into an output row" matrix product.
 */
export function axpy_simd_f64(out_ptr: usize, b_ptr: usize, scale: f64, len: i32): void {
    const len_simd = len - 1;
    const scale_v = f64x2.splat(scale);

    let j = 0;
    for (; j < len_simd; j += 2) {
        const offset = j << 3;
        const b_v = v128.load(b_ptr + offset);
        const out_v = v128.load(out_ptr + offset);
        v128.store(out_ptr + offset, f64x2.add(out_v, f64x2.mul(scale_v, b_v)));
    }

    for (; j < len; ++j) {
        const offset = j << 3;
        const prev = load<f64>(out_ptr + offset);
        store<f64>(out_ptr + offset, prev + scale * load<f64>(b_ptr + offset));
    }
}

/** Writes `len` zeroes at `ptr`. */
export function zero_f64(ptr: usize, len: i32): void {
    for (let i = 0; i < len; ++i) {
        store<f64>(ptr + (i << 3), 0.0);
    }
}

/** Sign of `x` as -1, 0 or 1, with zero its own sign. */
export function sign_f64(x: f64): i32 {
    return x > 0.0 ? 1 : x < 0.0 ? -1 : 0;
}

/**
 * The adaptive-gain rule shared by the `t-SNE` and `TriMap` optimisers.
 *
 * A gain shrinks by a fifth while the gradient keeps its direction and grows by 0.2 when the
 * direction flips, with a floor so a coordinate can never stop moving altogether.
 *
 * Zero counts as its own sign, so a zero paired with a non-zero is *not* the same sign. Both
 * callers start their step buffer at all zeros, and treating them as matching would take the wrong
 * branch across the whole first iteration.
 */
export function adapt_gain(gain: f64, grad: f64, step: f64): f64 {
    const same = sign_f64(grad) == sign_f64(step);
    const updated = same ? gain * 0.8 : gain + 0.2;
    return updated < 0.01 ? 0.01 : updated;
}
