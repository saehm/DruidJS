/**
 * WASM wrappers for the `SAMMON` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `SAMMON.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { copy_once, get_session, initWasm } from "../wasm/index.js";

/**
 * Sammon mapping single iteration gradient step using WASM SIMD.
 *
 * Buffers persist across calls, so the N ⨯ N distance matrix is copied in only on the first
 * iteration of a run. The caller must pass the same `D_val` array every iteration; see
 * {@link copy_once}. The gradient buffer is scratch the kernel fills before reading, so JS neither
 * provides nor receives it.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} D_val - Constant for the run; copied on first use only.
 * @param {number} n
 * @param {number} dim
 * @param {number} magic
 * @returns {boolean}
 */
export function wasmSammonStep(Y_val, D_val, n, dim, magic) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;

    const s = get_session(inst, exports, "sammon", {
        Y: Y_val.byteLength,
        D: D_val.byteLength,
        grad: n * dim * 8,
    });
    const { Y, D, grad } = s.ptrs;

    copy_once(exports, s, "D", D_val);

    const memory = exports.memory;
    new Float64Array(memory.buffer, Y, Y_val.length).set(Y_val);

    exports.sammon_step_f64(Y, D, grad, n, dim, magic);

    Y_val.set(new Float64Array(memory.buffer, Y, Y_val.length));
    return true;
}
