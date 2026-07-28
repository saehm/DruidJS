/**
 * Throwaway prototype for WASM.md step 2, question 2.
 *
 * Question 1 (`wasm_resident_probe.js`) showed the single-pair crossover collapses when operands are
 * already in linear memory. This asks whether that survives in situ, in the structure that motivated
 * the whole line of work: does HNSW query time actually drop, at unchanged recall?
 *
 * Residency is faked exactly as the plan specifies -- one hand-allocated block holding every element
 * back to back, and a pointer tagged onto each element object. `HNSW.js` is untouched; the resident
 * path is injected purely through the `metric` parameter.
 *
 * Three metrics per configuration:
 *   js       - the library `cosine`, i.e. what HNSW uses today at these dimensions
 *   copy     - `cosine` forced through WASM on every call (threshold ignored), today's boundary cost
 *   resident - the kernel called on two pointers, no copy
 *
 * **Each variant runs in its own process.** HNSW calls the metric through one call site; running
 * three different metric functions in one process turns that site polymorphic and penalises whichever
 * ran last, which is enough to invert the result. The parent spawns children and aggregates medians.
 *
 * All three compute the same angular distance, so recall must come out identical. It is reported per
 * run; a divergence means the comparison is invalid, not that WASM is better.
 *
 *     node benchmark/wasm_resident_hnsw_probe.js
 */

import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import { HNSW } from "../src/knn/HNSW.js";
import { WASM_BASE64 } from "../src/wasm/wasm_bytes.js";

const K = 10;
const N = 5000;
const N_QUERIES = 200;
const DIMS = [8, 16, 32, 64, 128];
const REPEATS = 5;

/**
 * Self-contained xorshift32, matching `knn_full_grid_benchmark.js` so inputs never depend on the
 * library under test.
 *
 * @param {number} seed
 * @returns {() => number}
 */
function xorshift32(seed) {
    let x = seed | 0 || 1;
    return () => {
        x ^= x << 13;
        x ^= x >>> 17;
        x ^= x << 5;
        return (x >>> 0) / 4294967296;
    };
}

/**
 * @param {number} n
 * @param {number} d
 * @param {() => number} random
 * @returns {number[][]}
 */
function make_points(n, d, random) {
    return Array.from({ length: n }, () => Array.from({ length: d }, () => random()));
}

/**
 * Clamped `acos`, copied from `metrics/cosine.js` so every variant returns the same quantity and
 * the WASM columns are not flattered by skipping a transcendental the JS column pays for.
 *
 * @param {number} similarity
 * @returns {number}
 */
function angle_from_similarity(similarity) {
    if (!Number.isFinite(similarity)) return Math.PI / 2;
    return Math.acos(Math.min(1, Math.max(-1, similarity)));
}

/** The plain JS cosine, lifted out of the library so no threshold check sits in the hot loop. */
function cosine_js(a, b) {
    const n = a.length;
    let sum = 0;
    let sum_a = 0;
    let sum_b = 0;
    for (let i = 0; i < n; ++i) {
        sum += a[i] * b[i];
        sum_a += a[i] * a[i];
        sum_b += b[i] * b[i];
    }
    if (sum_a === 0 || sum_b === 0) return Math.PI / 2;
    return angle_from_similarity(sum / Math.sqrt(sum_a * sum_b));
}

/**
 * Exact K nearest neighbours under `cosine_js`, the recall reference.
 *
 * @param {number[][]} points
 * @param {number[][]} queries
 * @returns {Set<number>[]}
 */
function exact_knn(points, queries) {
    return queries.map((q) => {
        const scored = points.map((p, i) => ({ i, d: cosine_js(p, q) }));
        scored.sort((x, y) => x.d - y.d);
        return new Set(scored.slice(0, K).map((s) => s.i));
    });
}

/**
 * @param {Set<number>[]} truth
 * @param {number[][]} got
 * @returns {number} Recall@K in [0, 1].
 */
function recall_at(truth, got) {
    let hits = 0;
    for (let i = 0; i < truth.length; ++i) {
        for (const idx of got[i]) if (truth[i].has(idx)) ++hits;
    }
    return hits / (truth.length * K);
}

/**
 * @param {number[]} xs
 * @returns {number}
 */
function median(xs) {
    const s = [...xs].sort((a, b) => a - b);
    const m = s.length >> 1;
    return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}

// ---------------------------------------------------------------------------------------------
// Child: one variant, one dimension, one process.
// ---------------------------------------------------------------------------------------------

/**
 * @param {"js" | "copy" | "resident"} variant
 * @param {number} d
 */
function child(variant, d) {
    const random = xorshift32(42);
    const points = make_points(N, d, random);
    const queries = make_points(N_QUERIES, d, random);
    const truth = exact_knn(points, queries);

    const bytes = Buffer.from(WASM_BASE64, "base64");
    const instance = new WebAssembly.Instance(new WebAssembly.Module(bytes), {
        env: {
            abort: (_m, _f, l, c) => {
                throw new Error(`abort ${l}:${c}`);
            },
        },
    });
    const exports = instance.exports;

    // Fake residency: one block holding every element, plus a scratch slot for the query.
    const base = exports.allocate((N + 1) * d * 8);
    const block = new Float64Array(exports.memory.buffer, base, (N + 1) * d);
    for (let i = 0; i < N; ++i) block.set(points[i], i * d);
    const query_ptr = base + N * d * 8;
    const query_view = new Float64Array(exports.memory.buffer, query_ptr, d);

    let pts = points;
    let qs = queries;
    let metric = cosine_js;

    // Control: the JS metric, but over the same tagged arrays the resident variant must use. Tagging
    // moves an array out of V8's fast elements mode, and HNSW touches elements outside the metric
    // too, so without this the resident column is charged for the fake-residency trick itself.
    if (variant === "js_tagged") {
        pts = points.map((p) => {
            const a = p.slice();
            a.__ptr = 0;
            return a;
        });
        qs = queries.map((q) => {
            const a = q.slice();
            a.__ptr = 0;
            return a;
        });
    }

    if (variant === "copy") {
        metric = (a, b) => {
            const pa = exports.allocate(d * 8);
            const pb = exports.allocate(d * 8);
            try {
                new Float64Array(exports.memory.buffer, pa, d).set(a);
                new Float64Array(exports.memory.buffer, pb, d).set(b);
                return angle_from_similarity(1.0 - exports.cosine_distance_f64(pa, pb, d));
            } finally {
                exports.free(pa);
                exports.free(pb);
            }
        };
    } else if (variant === "resident") {
        // Tagging a named property onto an array moves it out of V8's fast elements mode, so the
        // tagged copies are used only by this variant.
        pts = points.map((p, i) => {
            const a = p.slice();
            a.__ptr = base + i * d * 8;
            return a;
        });
        qs = queries.map((q) => {
            const a = q.slice();
            a.__ptr = query_ptr;
            return a;
        });
        metric = (a, b) => angle_from_similarity(1.0 - exports.cosine_distance_f64(a.__ptr, b.__ptr, d));
    }

    const resident = variant === "resident";
    let calls = 0;
    const counted = (a, b) => {
        ++calls;
        return metric(a, b);
    };

    const t0 = performance.now();
    const index = new HNSW(pts, { metric: counted, m: 16, ef_construction: 200, ef: 50, seed: 1212 });
    const t_build = performance.now() - t0;

    // Warm up before timing: a cold query path against a warm one is the mistake called out in
    // WASM.md's working notes.
    for (const q of qs.slice(0, 20)) {
        if (resident) query_view.set(q);
        index.search(q, K);
    }

    /** @type {number[]} */
    const times = [];
    let got = [];
    let calls_per_query = 0;
    for (let r = 0; r < REPEATS; ++r) {
        calls = 0;
        const t1 = performance.now();
        got = [];
        for (const q of qs) {
            if (resident) query_view.set(q);
            got.push(index.search(q, K).map((x) => x.index));
        }
        times.push((performance.now() - t1) / qs.length);
        calls_per_query = calls / qs.length;
    }

    // Per-call metric cost on this same data. The result is accumulated into a sink that is printed,
    // so the call cannot be eliminated as dead code -- that mistake made an earlier version of this
    // script report 10 ns for a d = 16 JS cosine.
    let sink = 0;
    const a0 = pts[1];
    const b0 = qs[0];
    if (resident) query_view.set(queries[0]);
    for (let i = 0; i < 50000; ++i) sink += metric(a0, b0);
    const t2 = performance.now();
    for (let i = 0; i < 50000; ++i) sink += metric(a0, b0);
    const ns_per_call = ((performance.now() - t2) / 50000) * 1e6;

    process.stdout.write(
        `${JSON.stringify({
            variant,
            d,
            t_build,
            t_query: median(times),
            calls: calls_per_query,
            ns_per_call,
            recall: recall_at(truth, got),
            sink_finite: Number.isFinite(sink),
        })}\n`,
    );
}

// ---------------------------------------------------------------------------------------------
// Parent: spawn one child per (variant, dimension).
// ---------------------------------------------------------------------------------------------

const argv = process.argv.slice(2);
const variant_arg = argv.find((a) => a.startsWith("--variant="));

if (variant_arg) {
    child(variant_arg.split("=")[1], Number(argv.find((a) => a.startsWith("--dim=")).split("=")[1]));
} else {
    const self = fileURLToPath(import.meta.url);
    console.log(`HNSW, N = ${N}, K = ${K}, ${N_QUERIES} queries, metric = cosine, median of ${REPEATS}`);
    console.log("each variant in a fresh process; 'metric' = share of query time inside the metric");
    console.log(
        `${"d".padStart(5)} | ${"query js".padStart(9)} | ${"js+tag".padStart(9)} | ${"query copy".padStart(10)} | ${"query res".padStart(9)} | ${"res/tag".padStart(7)} | ${"calls".padStart(5)} | ${"js ns".padStart(6)} | ${"res ns".padStart(6)} | ${"metric".padStart(6)} | ${"limit".padStart(6)} | recall js/copy/res`,
    );

    for (const d of DIMS) {
        /** @type {Record<string, any>} */
        const out = {};
        for (const variant of ["js", "js_tagged", "copy", "resident"]) {
            const stdout = execFileSync(process.execPath, [self, `--variant=${variant}`, `--dim=${d}`], {
                encoding: "utf8",
            });
            out[variant] = JSON.parse(stdout.trim().split("\n").pop());
        }

        // What fraction of a JS query is spent inside the metric, and so the best a free metric
        // could possibly do.
        const share = (out.js.calls * out.js.ns_per_call) / (out.js.t_query * 1e6);
        const limit = 1 - share;

        console.log(
            `${String(d).padStart(5)} | ${out.js.t_query.toFixed(3).padStart(7)}ms | ${out.js_tagged.t_query.toFixed(3).padStart(7)}ms | ${out.copy.t_query.toFixed(3).padStart(8)}ms | ${out.resident.t_query.toFixed(3).padStart(7)}ms | ${(out.resident.t_query / out.js_tagged.t_query).toFixed(2).padStart(6)}× | ${out.js.calls.toFixed(0).padStart(5)} | ${out.js.ns_per_call.toFixed(0).padStart(6)} | ${out.resident.ns_per_call.toFixed(0).padStart(6)} | ${(share * 100).toFixed(0).padStart(5)}% | ${limit.toFixed(2).padStart(5)}× | ${out.js.recall.toFixed(3)}/${out.copy.recall.toFixed(3)}/${out.resident.recall.toFixed(3)}`,
        );
    }
}
