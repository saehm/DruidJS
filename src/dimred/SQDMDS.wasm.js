/**
 * WASM wrappers for the `SQDMDS` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `SQDMDS.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { alloc, copy_once, free_all, get_session, initWasm } from "../wasm/index.js";

/**
 * Computes SQDMDS quartet gradients in WASM SIMD.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} X_val
 * @param {Uint32Array} quartets_val
 * @param {Float64Array} grads_val
 * @param {number} n
 * @param {number} d_ld
 * @param {number} d_hd
 * @param {boolean} use_exaggeration
 * @param {boolean} is_precomputed
 * @returns {boolean}
 */
export function wasmSqdmdsFillGrads(
    Y_val,
    X_val,
    quartets_val,
    grads_val,
    n,
    d_ld,
    d_hd,
    use_exaggeration,
    is_precomputed,
) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;

    const num_quartets = quartets_val.length / 4;

    const s = get_session(inst, exports, "sqdmds_fill_grads", {
        Y: Y_val.byteLength,
        X: X_val.byteLength,
        quartets: quartets_val.byteLength,
        grads: grads_val.byteLength,
    });
    const { Y, X, quartets, grads } = s.ptrs;

    // The high-dimensional data is fixed for the run; the quartets are resampled every iteration.
    copy_once(exports, s, "X", X_val);

    const memory = exports.memory;
    new Float64Array(memory.buffer, Y, Y_val.length).set(Y_val);
    new Int32Array(memory.buffer, quartets, quartets_val.length).set(quartets_val);

    exports.sqdmds_fill_grads_f64(Y, X, quartets, num_quartets, grads, n, d_ld, d_hd, use_exaggeration, is_precomputed);

    grads_val.set(new Float64Array(memory.buffer, grads, grads_val.length));
    return true;
}

/**
 * Range-based SQDMDS quartet gradient calculation for multi-threaded worker execution.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} X_val
 * @param {Uint32Array} quartets_val
 * @param {Float64Array} grads_val
 * @param {number} n
 * @param {number} d_ld
 * @param {number} d_hd
 * @param {boolean} use_exaggeration
 * @param {boolean} is_precomputed
 * @param {number} start_q
 * @param {number} end_q
 * @returns {boolean}
 */
export function wasmSqdmdsFillGradsRange(
    Y_val,
    X_val,
    quartets_val,
    grads_val,
    n,
    d_ld,
    d_hd,
    use_exaggeration,
    is_precomputed,
    start_q,
    end_q,
) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;
    const memory = exports.memory;

    const num_quartets = quartets_val.length / 4;

    /** @type {number[]} */
    const ptrs = [];
    try {
        const ptrY = alloc(exports, ptrs, Y_val.byteLength);
        const ptrX = alloc(exports, ptrs, X_val.byteLength);
        const ptrQ = alloc(exports, ptrs, quartets_val.byteLength);
        const ptrG = alloc(exports, ptrs, grads_val.byteLength);

        new Float64Array(memory.buffer, ptrY, Y_val.length).set(Y_val);
        new Float64Array(memory.buffer, ptrX, X_val.length).set(X_val);
        new Int32Array(memory.buffer, ptrQ, quartets_val.length).set(quartets_val);

        exports.sqdmds_fill_grads_range_f64(
            ptrY,
            ptrX,
            ptrQ,
            num_quartets,
            ptrG,
            n,
            d_ld,
            d_hd,
            use_exaggeration,
            is_precomputed,
            start_q,
            end_q,
        );

        grads_val.set(new Float64Array(memory.buffer, ptrG, grads_val.length));
        return true;
    } finally {
        free_all(exports, ptrs);
    }
}

/**
 * Applies SQDMDS Nesterov momentum update and updates embedding Y in WASM SIMD.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} M_val
 * @param {Float64Array} Grads_val
 * @param {number} n
 * @param {number} d
 * @param {number} lr
 * @returns {boolean}
 */
export function wasmSqdmdsNestrovStep(Y_val, M_val, Grads_val, n, d, lr) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;

    // Every operand here is per-iteration state and all of it is N ⨯ d, so nothing can be skipped;
    // the session only saves the allocate/free pair on each of the three buffers.
    const s = get_session(inst, exports, "sqdmds_nestrov", {
        Y: Y_val.byteLength,
        M: M_val.byteLength,
        grads: Grads_val.byteLength,
    });
    const { Y, M, grads } = s.ptrs;

    const memory = exports.memory;
    new Float64Array(memory.buffer, Y, Y_val.length).set(Y_val);
    new Float64Array(memory.buffer, M, M_val.length).set(M_val);
    new Float64Array(memory.buffer, grads, Grads_val.length).set(Grads_val);

    exports.sqdmds_nestrov_step_f64(Y, M, grads, n, d, lr);

    Y_val.set(new Float64Array(memory.buffer, Y, Y_val.length));
    M_val.set(new Float64Array(memory.buffer, M, M_val.length));
    return true;
}
