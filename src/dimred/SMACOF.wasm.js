/**
 * WASM wrappers for the `SMACOF` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `SMACOF.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { copy_once, get_session, initWasm } from "../wasm/index.js";

/**
 * Executes a single SMACOF Guttman Transform iteration step in WASM.
 *
 * @param {Float64Array} Z_val
 * @param {Float64Array} TargetD_val
 * @param {Float64Array} ZNew_val
 * @param {number} n
 * @param {number} d
 * @returns {number | null} Computed stress value, or null if WASM is unavailable.
 */
export function wasmSmacofStep(Z_val, TargetD_val, ZNew_val, n, d) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;

    const s = get_session(inst, exports, "smacof", {
        Z: Z_val.byteLength,
        TargetD: TargetD_val.byteLength,
        ZNew: ZNew_val.byteLength,
    });
    const { Z, TargetD, ZNew } = s.ptrs;

    // The target distances are fixed for the whole run; only the configuration changes per step.
    copy_once(exports, s, "TargetD", TargetD_val);

    const memory = exports.memory;
    new Float64Array(memory.buffer, Z, Z_val.length).set(Z_val);

    const stress = exports.smacof_step_f64(Z, TargetD, ZNew, n, d);

    ZNew_val.set(new Float64Array(memory.buffer, ZNew, ZNew_val.length));
    return stress;
}
