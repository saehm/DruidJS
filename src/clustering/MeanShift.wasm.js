/**
 * WASM wrappers for the `MeanShift` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `MeanShift.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, free_all, initWasm } from "../wasm/index.js";

/**
 * Computes a single MeanShift iteration step in WASM SIMD.
 *
 * @param {Float64Array} points_val
 * @param {Float64Array} out_points_val
 * @param {number} n
 * @param {number} d
 * @param {number} bandwidth
 * @param {boolean} use_gaussian
 * @returns {number | null} Max shift magnitude, or null if WASM is unavailable.
 */
export function wasmMeanShiftStep(points_val, out_points_val, n, d, bandwidth, use_gaussian) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrIn = alloc(exports, ptrs, points_val.byteLength);
        const ptrOut = alloc(exports, ptrs, out_points_val.byteLength);

        new Float64Array(memory.buffer, ptrIn, points_val.length).set(points_val);

        const max_shift = exports.meanshift_step_f64(ptrIn, ptrOut, n, d, bandwidth, use_gaussian);

        out_points_val.set(new Float64Array(memory.buffer, ptrOut, out_points_val.length));
        return max_shift;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Range-based MeanShift iteration step for multi-threaded worker execution.
 *
 * @param {Float64Array} points_val
 * @param {Float64Array} out_points_val
 * @param {number} n
 * @param {number} d
 * @param {number} bandwidth
 * @param {boolean} use_gaussian
 * @param {number} start_i
 * @param {number} end_i
 * @returns {number | null}
 */
export function wasmMeanShiftStepRange(points_val, out_points_val, n, d, bandwidth, use_gaussian, start_i, end_i) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrIn = alloc(exports, ptrs, points_val.byteLength);
        const ptrOut = alloc(exports, ptrs, out_points_val.byteLength);

        new Float64Array(memory.buffer, ptrIn, points_val.length).set(points_val);

        const max_shift = exports.meanshift_step_range_f64(
            ptrIn,
            ptrOut,
            n,
            d,
            bandwidth,
            use_gaussian,
            start_i,
            end_i,
        );

        // Slice only, so concurrent workers sharing `out_points_val` do not clobber each other.
        const offset = start_i * d;
        const count = (end_i - start_i) * d;
        out_points_val.set(new Float64Array(memory.buffer, ptrOut + offset * 8, count), offset);
        return max_shift;
    } finally {
        free_all(exports, ptrs);
    }
}
