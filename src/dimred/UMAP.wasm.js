/**
 * WASM wrappers for the `UMAP` kernels.
 *
 * Each one copies its operands into linear memory, calls the kernel compiled from
 * `UMAP.as.ts`, and copies the result back, reporting failure so the caller can fall back to
 * its JS implementation. The shared runtime -- instance setup, allocation and the persistent
 * buffer sessions -- lives in `src/wasm/index.js`.
 *
 * @module
 */

import { copy_once, get_session, initWasm } from "../wasm/index.js";

/**
 * Runs one UMAP SGD epoch: attractive forces over the active edges plus repulsive forces over their
 * negative samples.
 *
 * `negative_samples` holds the node indices drawn by the caller's seeded randomizer, in edge order —
 * the RNG stays in JS so both paths consume the identical sequence. `embedding`, `epoch_of_next_sample`
 * and `epoch_of_next_negative_sample` are updated in place.
 *
 * @param {Float64Array} embedding - The `n_points ⨯ dim` embedding, mutated in place.
 * @param {Int32Array} head - Edge source indices.
 * @param {Int32Array} tail - Edge target indices.
 * @param {Float32Array} epochs_per_sample
 * @param {Float32Array} epoch_of_next_sample - Mutated in place.
 * @param {Float32Array} epochs_per_negative_sample
 * @param {Float32Array} epoch_of_next_negative_sample - Mutated in place.
 * @param {Int32Array} negative_samples - Pre-drawn negative sample node indices.
 * @param {number} dim
 * @param {number} iter - Current epoch.
 * @param {number} a
 * @param {number} b
 * @param {number} gamma - Repulsion strength.
 * @param {number} alpha - Current learning rate.
 * @returns {boolean} True if successfully executed via WASM.
 */
export function wasmUmapOptimizeEpoch(
    embedding,
    head,
    tail,
    epochs_per_sample,
    epoch_of_next_sample,
    epochs_per_negative_sample,
    epoch_of_next_negative_sample,
    negative_samples,
    dim,
    iter,
    a,
    b,
    gamma,
    alpha,
) {
    const inst = initWasm();
    if (!inst) return false;

    /** @type {any} */
    const exports = inst.exports;

    // The edge list and the two per-sample schedules never change after `init()`, so they are copied
    // once per run rather than once per epoch — four edge-sized arrays out of the eight.
    // `epoch_of_next_*` cannot join them: the caller reads both every epoch to draw its negative
    // samples, so they have to come back out.
    const s = get_session(inst, exports, "umap_epoch", {
        Y: embedding.byteLength,
        head: head.byteLength,
        tail: tail.byteLength,
        eps: epochs_per_sample.byteLength,
        eons: epoch_of_next_sample.byteLength,
        epns: epochs_per_negative_sample.byteLength,
        eonns: epoch_of_next_negative_sample.byteLength,
    });
    const { Y, head: ptrHead, tail: ptrTail, eps, eons, epns, eonns } = s.ptrs;

    copy_once(exports, s, "head", head);
    copy_once(exports, s, "tail", tail);
    copy_once(exports, s, "eps", epochs_per_sample);
    copy_once(exports, s, "epns", epochs_per_negative_sample);

    // The negative sample count varies from epoch to epoch, so it gets its own slot: sharing one
    // would resize the session every epoch and evict the four constants above. Capacity is rounded
    // up to a power of two so the reallocation stops happening after the first few epochs.
    const neg_bytes = Math.max(negative_samples.byteLength, 4);
    const neg_capacity = 2 ** Math.ceil(Math.log2(neg_bytes));
    const neg_session = get_session(inst, exports, "umap_neg", { neg: neg_capacity });
    const ptrNeg = neg_session.ptrs.neg;

    const memory = exports.memory;
    new Float64Array(memory.buffer, Y, embedding.length).set(embedding);
    new Float32Array(memory.buffer, eons, epoch_of_next_sample.length).set(epoch_of_next_sample);
    new Float32Array(memory.buffer, eonns, epoch_of_next_negative_sample.length).set(epoch_of_next_negative_sample);
    if (negative_samples.length > 0) {
        new Int32Array(memory.buffer, ptrNeg, negative_samples.length).set(negative_samples);
    }

    exports.umap_optimize_epoch_f64(
        Y,
        ptrHead,
        ptrTail,
        eps,
        eons,
        epns,
        eonns,
        ptrNeg,
        head.length,
        negative_samples.length,
        dim,
        iter,
        a,
        b,
        gamma,
        alpha,
    );

    embedding.set(new Float64Array(memory.buffer, Y, embedding.length));
    epoch_of_next_sample.set(new Float32Array(memory.buffer, eons, epoch_of_next_sample.length));
    epoch_of_next_negative_sample.set(new Float32Array(memory.buffer, eonns, epoch_of_next_negative_sample.length));
    return true;
}
