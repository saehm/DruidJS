/**
 * WASM wrappers for the {@link StressMDS} kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from `StressMDS.as.ts`,
 * and copies the result back, reporting failure so the caller can fall back to its JS
 * implementation. The shared runtime -- instance setup, allocation and the persistent buffer
 * sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { copy_once, get_session, initWasm } from "../wasm/index.js";

/**
 * Session buffer sizes. Both kernels share one session, so the shapes must agree between them.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} D_val
 * @returns {Record<string, number>}
 */
function sizes(Y_val, D_val) {
    return {
        Y: Y_val.byteLength,
        D: D_val.byteLength,
        W: D_val.byteLength,
        G: Y_val.byteLength,
        Trial: Y_val.byteLength,
        Out: 16,
    };
}

/**
 * Runs one weighted-stress iteration -- preconditioned gradient, line search and update -- in WASM.
 *
 * `Y_val` is updated in place when the step is accepted. The distances, weights, gradient and trial
 * scratch stay in the session between iterations, so a run copies the embedding in and out once per
 * iteration rather than once per line-search trial.
 *
 * @param {Float64Array} Y_val - N ⨯ d embedding, updated in place on an accepted step.
 * @param {Float64Array} D_val - N ⨯ N target distances. Constant for a run.
 * @param {Float64Array} W_val - N ⨯ N weights. Constant for a run.
 * @param {number} n
 * @param {number} d
 * @param {number} step - Current step size.
 * @param {number} energy_current - Energy of `Y_val` on entry.
 * @returns {{ energy: number; step: number; accepted: boolean } | null} The new state, or null if
 *   WASM is unavailable.
 */
export function wasmStressMDSStep(Y_val, D_val, W_val, n, d, step, energy_current) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;

    const s = get_session(inst, exports, "stress_mds", sizes(Y_val, D_val));
    const { Y, D, W, G, Trial, Out } = s.ptrs;

    // Distances and weights are fixed for the whole run; only the configuration changes per step.
    copy_once(exports, s, "D", D_val);
    copy_once(exports, s, "W", W_val);

    const memory = exports.memory;
    new Float64Array(memory.buffer, Y, Y_val.length).set(Y_val);

    const energy = exports.stress_mds_step_f64(Y, D, W, G, Trial, n, d, step, energy_current, Out);

    const out = new Float64Array(memory.buffer, Out, 2);
    const accepted = out[1] === 1;
    if (accepted) {
        Y_val.set(new Float64Array(memory.buffer, Y, Y_val.length));
    }
    return { energy, step: out[0], accepted };
}

/**
 * Computes the weighted stress of a configuration in WASM.
 *
 * Called once to seed the iteration loop; later values come back from {@link wasmStressMDSStep}.
 *
 * @param {Float64Array} Y_val
 * @param {Float64Array} D_val
 * @param {Float64Array} W_val
 * @param {number} n
 * @param {number} d
 * @returns {number | null} The stress, or null if WASM is unavailable.
 */
export function wasmStressMDSEnergy(Y_val, D_val, W_val, n, d) {
    const inst = initWasm();
    if (!inst) return null;

    /** @type {any} */
    const exports = inst.exports;

    const s = get_session(inst, exports, "stress_mds", sizes(Y_val, D_val));
    const { Y, D, W } = s.ptrs;

    copy_once(exports, s, "D", D_val);
    copy_once(exports, s, "W", W_val);
    new Float64Array(exports.memory.buffer, Y, Y_val.length).set(Y_val);

    return exports.stress_mds_energy_f64(Y, D, W, n, d);
}
