/**
 * WASM wrappers for the `TSNE` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `TSNE.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { copy_once, get_session, initWasm } from "../wasm/index.js";

/**
 * Performs a single t-SNE iteration gradient and position update step via WASM.
 *
 * Buffers persist across calls, so `P` is copied in only on the first iteration of a run and the
 * N ⨯ N scratch buffers — which JS neither writes nor reads — are never copied at all. The caller
 * must therefore pass the same `P_val` array every iteration; a different array is detected and
 * recopied, but mutating one in place behind a stable identity is not. See {@link copy_once}.
 *
 * `Y`, `ystep` and `gains` are N ⨯ dim, small next to the N ⨯ N operands, and are still copied both
 * ways each call so that interleaved runs of the same size cannot read each other's state.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} P_val - Constant for the run; copied on first use only.
 * @param {Float64Array} ystep_val
 * @param {Float64Array} gains_val
 * @param {number} n
 * @param {number} dim
 * @param {number} pmul
 * @param {number} epsilon
 * @param {number} momval
 * @returns {boolean} True if executed successfully via WASM.
 */
export function wasmTsneStep(Y_val, P_val, ystep_val, gains_val, n, dim, pmul, epsilon, momval) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;

    // Buffers backed by a JS array are sized from that array; only the pure scratch, which JS never
    // sees, is sized from the kernel's own indexing.
    const s = get_session(inst, exports, "tsne", {
        Y: Y_val.byteLength,
        P: P_val.byteLength,
        Qu: n * n * 8,
        Q: n * n * 8,
        grad: n * dim * 8,
        ystep: ystep_val.byteLength,
        gains: gains_val.byteLength,
    });
    const { Y, P, Qu, Q, grad, ystep, gains } = s.ptrs;

    copy_once(exports, s, "P", P_val);

    // `memory.buffer` is re-read after the copies above: allocating can grow linear memory, which
    // replaces the buffer and detaches every view onto the old one.
    const memory = exports.memory;
    new Float64Array(memory.buffer, Y, Y_val.length).set(Y_val);
    new Float64Array(memory.buffer, ystep, ystep_val.length).set(ystep_val);
    new Float64Array(memory.buffer, gains, gains_val.length).set(gains_val);

    exports.tsne_step_f64(Y, P, Qu, Q, grad, ystep, gains, n, dim, pmul, epsilon, momval);

    Y_val.set(new Float64Array(memory.buffer, Y, Y_val.length));
    ystep_val.set(new Float64Array(memory.buffer, ystep, ystep_val.length));
    gains_val.set(new Float64Array(memory.buffer, gains, gains_val.length));
    return true;
}
