import { wasmCosineDistance } from "../wasm/index.js";
import { WASM_MIN_VECTOR_LENGTH } from "../wasm/thresholds.js";

/**
 * Turns a cosine similarity into an angular distance.
 *
 * The similarity is clamped to `[-1, 1]` before taking the arc cosine: floating point rounding
 * lets the ratio drift marginally outside that range (e.g. `1.0000000000000002` for two identical
 * vectors), which would otherwise make `Math.acos` return `NaN`.
 *
 * A vector of length zero has no direction, so the angle to it is undefined. Rather than
 * propagating `NaN` through a whole embedding, both the JS and the WASM path report such pairs as
 * orthogonal (`π / 2`).
 *
 * @private
 * @param {number} similarity - The cosine similarity, nominally in `[-1, 1]`.
 * @returns {number} The angle in radians.
 */
function angle_from_similarity(similarity) {
    if (!Number.isFinite(similarity)) return Math.PI / 2;
    return Math.acos(Math.min(1, Math.max(-1, similarity)));
}

/**
 * Computes the cosine distance (not similarity) between `a` and `b`.
 *
 * Identical vectors yield exactly `0`, and pairs involving a zero vector yield `π / 2`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The cosine distance between `a` and `b`.
 * @example
 * import { cosine } from "@saehrimnir/druidjs";
 * const a = [1, 2, 3];
 * const b = [4, 5, 6];
 * const distance = cosine(a, b); // 0.9746318461970762
 */
export function cosine(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;

    if (n >= WASM_MIN_VECTOR_LENGTH) {
        const wasmRes = wasmCosineDistance(a, b);
        if (wasmRes !== null) {
            return angle_from_similarity(1.0 - wasmRes);
        }
    }

    let sum = 0;
    let sum_a = 0;
    let sum_b = 0;
    for (let i = 0; i < n; ++i) {
        sum += a[i] * b[i];
        sum_a += a[i] * a[i];
        sum_b += b[i] * b[i];
    }
    if (sum_a === 0 || sum_b === 0) return Math.PI / 2;
    // One square root of the product, not the product of two square roots: `Math.sqrt(s) *
    // Math.sqrt(s)` rounds twice and lands either side of `s`, so identical vectors would yield a
    // similarity a few ulps off 1 and an angle of ~1e-8 instead of 0.
    return angle_from_similarity(sum / Math.sqrt(sum_a * sum_b));
}
