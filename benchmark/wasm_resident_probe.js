/**
 * Throwaway prototype for WASM.md step 2, question 1.
 *
 * Where is the break-even for a single-pair metric call when the operands are ALREADY in WASM
 * linear memory? Residency is faked with one hand-allocated block, as the plan specifies -- no
 * Matrix change. If this does not push the crossover below ~64 the step-2 premise is dead.
 *
 * Three variants per size:
 *   js       - the JS loop over JS-heap Float64Arrays (what the library does today below threshold)
 *   copy     - the current WASM path: allocate, copy both vectors in, call, free
 *   resident - the step-2 path: both vectors already in linear memory, pass pointers only
 *
 * Pairs are drawn from a pool of many vectors rather than reusing one pair, so the numbers reflect
 * KNN-style access (cache-cold-ish, no JIT hoisting) instead of a hot two-vector loop.
 */

import { WASM_BASE64 } from "../src/wasm/wasm_bytes.js";

const bytes = Buffer.from(WASM_BASE64, "base64");
const instance = new WebAssembly.Instance(new WebAssembly.Module(bytes), {
    env: { abort: (_m, _f, l, c) => { throw new Error(`abort ${l}:${c}`); } },
});
const exports = instance.exports;

/** Number of vectors in the pool; enough that the working set exceeds L2 at larger d. */
const POOL = 256;

/**
 * @param {number} n
 * @returns {Float64Array}
 */
function random_vector(n) {
    const a = new Float64Array(n);
    for (let i = 0; i < n; ++i) a[i] = Math.random();
    return a;
}

/**
 * @param {() => number} fn
 * @param {number} reps
 * @returns {number} Mean nanoseconds per call.
 */
function bench(fn, reps) {
    let sink = 0;
    for (let i = 0; i < reps; ++i) sink += fn(); // warm up the JIT before timing
    const start = performance.now();
    for (let i = 0; i < reps; ++i) sink += fn();
    const ns = ((performance.now() - start) / reps) * 1e6;
    if (!Number.isFinite(sink)) throw new Error("sink went non-finite");
    return ns;
}

const SIZES = [4, 8, 16, 32, 64, 128, 256, 512, 1024];

// The floor under every resident call: crossing into WASM and back with nothing to compute.
// `len = 0` makes the kernel return on its zero-norm guard before touching memory, so this is the
// boundary alone. Residency removes the *copy*; it cannot remove this.
{
    const p = exports.allocate(8);
    let sink = 0;
    for (let i = 0; i < 200000; ++i) sink += exports.cosine_distance_f64(p, p, 0);
    const t = performance.now();
    for (let i = 0; i < 200000; ++i) sink += exports.cosine_distance_f64(p, p, 0);
    const ns = ((performance.now() - t) / 200000) * 1e6;
    exports.free(p);
    if (!Number.isFinite(sink)) throw new Error("sink went non-finite");
    console.log(`empty JS->WASM call (len = 0): ${ns.toFixed(1)} ns -- the floor for any resident call`);
}

/**
 * @param {string} title
 * @param {string} kernel - Exported kernel name.
 * @param {(a: Float64Array, b: Float64Array, n: number) => number} js_loop
 */
function probe(title, kernel, js_loop, post) {
    console.log(`\n${title}`);
    console.log(
        `${"d".padStart(6)} | ${"js (ns)".padStart(9)} | ${"copy (ns)".padStart(9)} | ${"resident (ns)".padStart(13)} | ${"tagged".padStart(10)} | ${"copy/js".padStart(8)} | resident/js`,
    );
    for (const n of SIZES) {
        // JS-heap pool.
        const heap_pool = [];
        for (let i = 0; i < POOL; ++i) heap_pool.push(random_vector(n));

        // Resident pool: one allocation holding every vector back to back.
        const base = exports.allocate(POOL * n * 8);
        const view = new Float64Array(exports.memory.buffer, base, POOL * n);
        for (let i = 0; i < POOL; ++i) view.set(heap_pool[i], i * n);

        // Deterministic, non-sequential pair walk. Each variant keeps its own counter and writes the
        // two indices to module-level slots -- no allocation, so the harness does not bias small d.
        const reps = Math.max(2000, Math.min(50000, (1 << 20) / n));

        // What the walk alone costs, so it can be subtracted from every column.
        let k_null = 0;
        const t_null = bench(() => {
            k_null = (k_null + 1) % (POOL * POOL);
            return ((k_null * 37) % POOL) + ((k_null * 91 + 13) % POOL);
        }, reps);

        let k_js = 0;
        const t_js = bench(() => {
            k_js = (k_js + 1) % (POOL * POOL);
            const i = (k_js * 37) % POOL;
            const j = (k_js * 91 + 13) % POOL;
            return js_loop(heap_pool[i], heap_pool[j], n);
        }, reps);

        let k_copy = 0;
        const t_copy = bench(() => {
            k_copy = (k_copy + 1) % (POOL * POOL);
            const i = (k_copy * 37) % POOL;
            const j = (k_copy * 91 + 13) % POOL;
            const pa = exports.allocate(n * 8);
            const pb = exports.allocate(n * 8);
            try {
                const ma = new Float64Array(exports.memory.buffer, pa, n);
                const mb = new Float64Array(exports.memory.buffer, pb, n);
                ma.set(heap_pool[i]);
                mb.set(heap_pool[j]);
                return post(exports[kernel](pa, pb, n));
            } finally {
                exports.free(pa);
                exports.free(pb);
            }
        }, reps);

        let k_res = 0;
        const t_res = bench(() => {
            k_res = (k_res + 1) % (POOL * POOL);
            const i = (k_res * 37) % POOL;
            const j = (k_res * 91 + 13) % POOL;
            return post(exports[kernel](base + i * n * 8, base + j * n * 8, n));
        }, reps);

        // Same call, but the pointer arrives as a named property on a tagged array -- how the HNSW
        // probe has to fake residency. Isolates what the tag costs from what residency buys.
        const tagged = [];
        for (let i = 0; i < POOL; ++i) {
            const a = Array.from(heap_pool[i]);
            a.__ptr = base + i * n * 8;
            tagged.push(a);
        }
        let k_tag = 0;
        const t_tag = bench(() => {
            k_tag = (k_tag + 1) % (POOL * POOL);
            const i = (k_tag * 37) % POOL;
            const j = (k_tag * 91 + 13) % POOL;
            return post(exports[kernel](tagged[i].__ptr, tagged[j].__ptr, n));
        }, reps);

        exports.free(base);

        // Subtract the walk so the columns are the metric call itself.
        const js = t_js - t_null;
        const copy = t_copy - t_null;
        const res = t_res - t_null;

        const tag = t_tag - t_null;

        const mark = res < js ? "" : "  (JS wins)";
        console.log(
            `${String(n).padStart(6)} | ${js.toFixed(1).padStart(9)} | ${copy.toFixed(1).padStart(9)} | ${res.toFixed(1).padStart(13)} | ${tag.toFixed(1).padStart(10)} | ${(copy / js).toFixed(2).padStart(7)}× | ${(res / js).toFixed(2)}×${mark}`,
        );
    }
}

probe(
    "cosine — multi-accumulator (WASM_MIN_VECTOR_LENGTH = 512)",
    "cosine_distance_f64",
    (a, b, n) => {
        let sum = 0;
        let sum_a = 0;
        let sum_b = 0;
        for (let i = 0; i < n; ++i) {
            sum += a[i] * b[i];
            sum_a += a[i] * a[i];
            sum_b += b[i] * b[i];
        }
        return Math.acos(Math.min(1, Math.max(-1, sum / Math.sqrt(sum_a * sum_b))));
    },
    // `metrics/cosine.js` turns the kernel's `1 - similarity` into an angle. Charging that to the JS
    // column but not the WASM ones overstated the win by roughly the cost of an `acos`.
    (raw) => Math.acos(Math.min(1, Math.max(-1, 1.0 - raw))),
);

probe(
    "manhattan — single-accumulator (WASM_MIN_SIMPLE_VECTOR_LENGTH = 1024)",
    "manhattan_distance_f64",
    (a, b, n) => {
        let sum = 0;
        for (let i = 0; i < n; ++i) sum += Math.abs(a[i] - b[i]);
        return sum;
    },
    // `metrics/manhattan.js` returns the kernel result unchanged.
    (raw) => raw,
);
