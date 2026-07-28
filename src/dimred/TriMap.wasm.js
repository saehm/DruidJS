/**
 * WASM wrappers for the `TriMap` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `TriMap.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { copy_once, get_session, initWasm } from "../wasm/index.js";

/**
 * Computes TriMap triplet gradients in WASM SIMD.
 *
 * Buffers persist across calls, so the triplet list and their weights are copied in only once per
 * run. They are the largest operands here — there are far more triplets than points — and `init()`
 * builds both once, so the caller must pass the same arrays every iteration; see {@link copy_once}.
 * `Y` and `grad` are rebuilt by the caller each iteration and are copied as before.
 *
 * @param {Float64Array} Y_val
 * @param {Int32Array} triplets_val - Constant for the run; copied on first use only.
 * @param {Float64Array} weights_val - Constant for the run; copied on first use only.
 * @param {Float64Array} grad_val
 * @param {number} n
 * @param {number} dim
 * @param {number} n_inliers
 * @param {number} n_outliers
 * @returns {number | null} Computed loss, or null if WASM is unavailable.
 */
export function wasmTriMapGrad(Y_val, triplets_val, weights_val, grad_val, n, dim, n_inliers, n_outliers) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;

    const num_triplets = triplets_val.length / 3;

    const s = get_session(inst, exports, "trimap_grad", {
        Y: Y_val.byteLength,
        triplets: triplets_val.byteLength,
        weights: weights_val.byteLength,
        grad: grad_val.byteLength,
    });
    const { Y, triplets, weights, grad } = s.ptrs;

    copy_once(exports, s, "triplets", triplets_val);
    copy_once(exports, s, "weights", weights_val);

    const memory = exports.memory;
    new Float64Array(memory.buffer, Y, Y_val.length).set(Y_val);

    const loss = exports.trimap_grad_f64(Y, triplets, weights, grad, n, dim, num_triplets, n_inliers, n_outliers);

    grad_val.set(new Float64Array(memory.buffer, grad, grad_val.length));
    return loss;
}

/**
 * Applies TriMap embedding update in WASM SIMD.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} grad_val
 * @param {Float64Array} vel_val
 * @param {Float64Array} gain_val
 * @param {number} n
 * @param {number} dim
 * @param {number} gamma
 * @param {number} lr
 * @returns {boolean}
 */
export function wasmTriMapUpdate(Y_val, grad_val, vel_val, gain_val, n, dim, gamma, lr) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;

    // All four operands are per-iteration state and all are N ⨯ dim, so nothing can be skipped here;
    // the session only saves the allocate/free pair on each buffer.
    const s = get_session(inst, exports, "trimap_update", {
        Y: Y_val.byteLength,
        grad: grad_val.byteLength,
        vel: vel_val.byteLength,
        gain: gain_val.byteLength,
    });
    const { Y, grad, vel, gain } = s.ptrs;

    const memory = exports.memory;
    new Float64Array(memory.buffer, Y, Y_val.length).set(Y_val);
    new Float64Array(memory.buffer, grad, grad_val.length).set(grad_val);
    new Float64Array(memory.buffer, vel, vel_val.length).set(vel_val);
    new Float64Array(memory.buffer, gain, gain_val.length).set(gain_val);

    exports.trimap_update_f64(Y, grad, vel, gain, n, dim, gamma, lr);

    Y_val.set(new Float64Array(memory.buffer, Y, Y_val.length));
    vel_val.set(new Float64Array(memory.buffer, vel, vel_val.length));
    gain_val.set(new Float64Array(memory.buffer, gain, gain_val.length));
    return true;
}
