import { describe, expect, it } from "vitest";
import { TSNE } from "../../src/dimred/TSNE.js";
import { distance_matrix } from "../../src/matrix/distance_matrix.js";
import { Matrix } from "../../src/matrix/Matrix.js";
import { held_session_count, isWasmAvailable, release_wasm_buffers, setWasmEnabled } from "../../src/wasm/index.js";

/**
 * The iterative kernels keep their WASM buffers alive between calls instead of allocating and
 * copying per iteration. That makes a whole class of mistake possible that the per-call code could
 * not have: a buffer reused at the wrong size, state surviving from a previous run, or a pointer
 * outliving the heap it came from.
 *
 * The first of those already happened. t-SNE's `ystep` and `gains` are N ⨯ D in the *input*
 * dimensionality, not the N ⨯ dim the kernel indexes, so sizing their buffers from `n * dim` while
 * viewing them at `ystep.length` wrote through the next allocator block header. Nothing failed at
 * the time — the damage only surfaced as an abort inside `~lib/rt/tlsf.ts` on the *next* unrelated
 * allocation, which is why "t-SNE then something else" is a test here rather than an assumption.
 */

/**
 * @param {number} n
 * @param {number} d
 * @returns {number[][]} Seeded so a failure is reproducible.
 */
function data(n, d) {
    let x = 42;
    const random = () => {
        x ^= x << 13;
        x ^= x >>> 17;
        x ^= x << 5;
        return (x >>> 0) / 4294967296;
    };
    return Array.from({ length: n }, () => Array.from({ length: d }, () => random()));
}

/**
 * @param {number[][]} X
 * @param {number} iterations
 * @returns {number[]}
 */
function project(X, iterations = 5) {
    return new TSNE(X, { d: 2, seed: 42, perplexity: 10 }).transform(iterations).flat();
}

/**
 * @param {number[]} values
 * @returns {boolean}
 */
function all_finite(values) {
    return values.length > 0 && values.every((v) => Number.isFinite(v));
}

describe.runIf(isWasmAvailable())("persistent WASM buffers", () => {
    it("leaves the heap intact for the next kernel", () => {
        const X = data(40, 5);
        project(X);

        // Allocates in the same heap. Before the size fix this aborted in the TLSF allocator.
        const D = distance_matrix(Matrix.from(X));

        expect(D.shape).toEqual([40, 40]);
        expect(Number.isFinite(D.entry(0, 1))).toBe(true);
        expect(isWasmAvailable(), "a trap would have retired the instance").toBe(true);
    });

    it("reallocates when the input dimensionality changes but n and dim do not", () => {
        // Same 40 ⨯ 2 output and so the same `n` and `dim`, but `ystep`/`gains` are n ⨯ D and D
        // differs. A session keyed on n and dim alone would reuse buffers that are too small here.
        const wide = project(data(40, 12));
        const narrow = project(data(40, 3));

        expect(all_finite(wide)).toBe(true);
        expect(all_finite(narrow)).toBe(true);
        expect(isWasmAvailable()).toBe(true);
    });

    it("does not carry state between runs of the same shape", () => {
        const X = data(40, 5);
        const first = project(X);
        const second = project(X);

        // Same input and seed: reusing the buffers must not make the second run differ from the
        // first, which it would if scratch or `P` leaked across runs.
        expect(second).toEqual(first);
    });

    it("interleaves runs of different sizes without corruption", () => {
        const small = data(20, 4);
        const large = data(60, 7);

        const small_alone = project(small);
        const large_alone = project(large);

        // Force the session to be torn down and rebuilt repeatedly.
        project(small);
        project(large);
        const small_again = project(small);
        const large_again = project(large);

        expect(small_again).toEqual(small_alone);
        expect(large_again).toEqual(large_alone);
    });

    it("keeps working after the buffers are released", () => {
        const X = data(40, 5);
        const before = project(X);

        release_wasm_buffers();
        release_wasm_buffers(); // idempotent — a second release must not double free

        expect(project(X)).toEqual(before);
    });

    it("releases its buffers when a run ends", () => {
        const X = data(40, 5);
        release_wasm_buffers();

        project(X);
        expect(held_session_count(), "transform() must not leave buffers behind").toBe(0);
    });

    it("releases its buffers when a generator completes", () => {
        const X = data(40, 5);
        release_wasm_buffers();

        for (const _ of new TSNE(X, { d: 2, seed: 42, perplexity: 10 }).generator(5)) {
            // drained
        }
        expect(held_session_count()).toBe(0);
    });

    it("releases its buffers when a generator is abandoned by `break`", () => {
        const X = data(40, 5);
        release_wasm_buffers();

        for (const _ of new TSNE(X, { d: 2, seed: 42, perplexity: 10 }).generator(50)) {
            break; // `for…of` calls return(), which must run the finally
        }
        expect(held_session_count()).toBe(0);
    });

    it("releases through release() for a hand-driven run", () => {
        const X = data(40, 5);
        release_wasm_buffers();

        const tsne = new TSNE(X, { d: 2, seed: 42, perplexity: 10 });
        const steps = tsne.generator(50);
        steps.next();
        steps.next();
        expect(held_session_count(), "a half-driven generator still holds its buffers").toBe(1);

        expect(tsne.release(), "release() returns the instance for chaining").toBe(tsne);
        expect(held_session_count()).toBe(0);
        tsne.release(); // idempotent — a second call must not throw or double free
        expect(held_session_count()).toBe(0);
    });

    it("release() on one instance leaves another run of the same method correct", () => {
        const X = data(40, 5);
        const reference = project(X);

        // Two hand-driven runs of the same method share the session slot (keyed by kernel, not
        // instance). Releasing the first must not corrupt the second — it only forces a realloc.
        const a = new TSNE(X, { d: 2, seed: 42, perplexity: 10 });
        const b = new TSNE(X, { d: 2, seed: 42, perplexity: 10 });
        const ga = a.generator(5);
        const gb = b.generator(5);
        ga.next();
        gb.next();
        a.release(); // frees the slot b is mid-run using

        // b finishes; its next step reallocates and re-copies its JS-held state, so the result must
        // still match a clean, uninterrupted run.
        let last = /** @type {number[][]} */ ([]);
        for (const step of gb) last = /** @type {number[][]} */ (step);
        expect(last.flat()).toEqual(reference);
    });

    it("releases through Symbol.dispose for a hand-driven run", () => {
        const X = data(40, 5);
        release_wasm_buffers();

        const tsne = new TSNE(X, { d: 2, seed: 42, perplexity: 10 });
        const steps = tsne.generator(50);
        steps.next();
        steps.next();
        expect(held_session_count(), "a half-driven generator still holds its buffers").toBe(1);

        tsne[Symbol.dispose]();
        expect(held_session_count()).toBe(0);
    });

    it("agrees with the JS fallback", () => {
        const X = data(40, 5);

        setWasmEnabled(true);
        const wasm = project(X, 20);
        setWasmEnabled(false);
        const js = project(X, 20);
        setWasmEnabled(true);

        // t-SNE is one of the paths that comes out bit-identical, so this is an exact comparison.
        expect(wasm).toEqual(js);
    });
});
