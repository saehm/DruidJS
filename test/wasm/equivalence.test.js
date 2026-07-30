import { afterAll, describe, expect, it } from "vitest";
import { KMeans } from "../../src/clustering/KMeans.js";
import { MeanShift } from "../../src/clustering/MeanShift.js";
import { ISOMAP } from "../../src/dimred/ISOMAP.js";
import { SAMMON } from "../../src/dimred/SAMMON.js";
import { SMACOF } from "../../src/dimred/SMACOF.js";
import { SQDMDS } from "../../src/dimred/SQDMDS.js";
import { TriMap } from "../../src/dimred/TriMap.js";
import { TSNE } from "../../src/dimred/TSNE.js";
import { wasmUmapOptimizeEpoch } from "../../src/dimred/UMAP.wasm.js";
import { inner_product } from "../../src/linear_algebra/inner_product.js";
import { distance_matrix } from "../../src/matrix/distance_matrix.js";
import { Matrix } from "../../src/matrix/Matrix.js";
import { norm } from "../../src/matrix/norm.js";
import { normalize } from "../../src/matrix/normalize.js";
import { bray_curtis } from "../../src/metrics/bray_curtis.js";
import { canberra } from "../../src/metrics/canberra.js";
import { chebyshev } from "../../src/metrics/chebyshev.js";
import { cosine } from "../../src/metrics/cosine.js";
import { manhattan } from "../../src/metrics/manhattan.js";
import { Randomizer } from "../../src/util/randomizer.js";
import { isWasmAvailable, setWasmEnabled } from "../../src/wasm/index.js";
import { WASM_MIN_UMAP_EDGES } from "../../src/wasm/thresholds.js";

/**
 * Every accelerated function has two implementations. The WASM kernel is the one that runs in CI on
 * a capable runner, which leaves the JS fallback — the code that actually runs in restricted
 * browsers — completely unexercised. These tests run both and assert they agree.
 *
 * The two are not bit-identical by construction: the kernels accumulate in SIMD lanes, so partial
 * sums are combined in a different order than the scalar loops. Agreement is asserted to a relative
 * tolerance rather than exactly.
 */

const RELATIVE_TOLERANCE = 1e-9;

/**
 * Runs `fn` twice, once with WASM enabled and once forced onto the JS path.
 *
 * @template T
 * @param {() => T} fn
 * @returns {{ wasm: T; js: T }}
 */
function both_paths(fn) {
    setWasmEnabled(true);
    const wasm = fn();
    setWasmEnabled(false);
    const js = fn();
    setWasmEnabled(true);
    return { wasm, js };
}

/**
 * @param {number} n
 * @param {number} d
 * @param {number} [seed]
 * @returns {number[][]}
 */
function random_data(n, d, seed = 1212) {
    const R = new Randomizer(seed);
    return Array.from({ length: n }, () => Array.from({ length: d }, () => R.random));
}

/**
 * @param {unknown} a
 * @param {unknown} b
 * @param {string} label
 */
function expect_close(a, b, label) {
    const flat_a = /** @type {number[]} */ ([]).concat(...[a].flat(3));
    const flat_b = /** @type {number[]} */ ([]).concat(...[b].flat(3));
    expect(flat_a.length, `${label}: length`).toBe(flat_b.length);
    for (let i = 0; i < flat_a.length; ++i) {
        const scale = Math.max(1, Math.abs(flat_a[i]), Math.abs(flat_b[i]));
        expect(Math.abs(flat_a[i] - flat_b[i]), `${label}: index ${i}`).toBeLessThan(RELATIVE_TOLERANCE * scale);
    }
}

afterAll(() => {
    setWasmEnabled(true);
});

describe.runIf(isWasmAvailable())("WASM and JS fallback agree", () => {
    it("exposes a working opt-out switch", () => {
        expect(isWasmAvailable()).toBe(true);
        setWasmEnabled(false);
        expect(isWasmAvailable()).toBe(false);
        setWasmEnabled(true);
        expect(isWasmAvailable()).toBe(true);
    });

    // Long enough to clear every vector threshold in src/wasm/thresholds.js.
    const a = Float64Array.from({ length: 2048 }, (_, i) => Math.sin(i) + 2);
    const b = Float64Array.from({ length: 2048 }, (_, i) => Math.cos(i) + 2);

    for (const [name, metric] of [
        ["cosine", cosine],
        ["manhattan", manhattan],
        ["chebyshev", chebyshev],
        ["canberra", canberra],
        ["bray_curtis", bray_curtis],
    ]) {
        it(`${name} agrees on both paths`, () => {
            const { wasm, js } = both_paths(() => metric(a, b));
            expect_close(wasm, js, name);
        });
    }

    it("norm agrees on both paths", () => {
        const { wasm, js } = both_paths(() => norm(a));
        expect_close(wasm, js, "norm");
    });

    it("normalize agrees on both paths", () => {
        const { wasm, js } = both_paths(() => Array.from(normalize(a)));
        expect_close(wasm, js, "normalize");
    });

    it("inner_product agrees on both paths", () => {
        const { wasm, js } = both_paths(() => inner_product(a, b));
        expect_close(wasm, js, "inner_product");
    });

    it("Matrix.dot / transDot / dotTrans / outer agree on both paths", () => {
        const A = new Matrix(24, 17, (i, j) => Math.sin(i * 17 + j));
        const B = new Matrix(17, 21, (i, j) => Math.cos(i * 21 + j));
        expect_close(...Object.values(both_paths(() => A.dot(B).asArray())), "dot");
        expect_close(...Object.values(both_paths(() => A.T.transDot(B).asArray())), "transDot");
        expect_close(...Object.values(both_paths(() => A.dotTrans(B.T).asArray())), "dotTrans");

        const v = Matrix.from_vector(
            Float64Array.from({ length: 40 }, (_, i) => Math.sin(i)),
            "col",
        );
        expect_close(...Object.values(both_paths(() => v.outer(v).asArray())), "outer");
    });

    it("distance_matrix agrees on both paths", () => {
        const X = random_data(40, 6);
        const { wasm, js } = both_paths(() => distance_matrix(Matrix.from(X)).asArray());
        expect_close(wasm, js, "distance_matrix");
    });

    it("distance_matrix stays symmetric with a zero diagonal", () => {
        const D = distance_matrix(Matrix.from(random_data(40, 6)));
        for (let i = 0; i < 40; ++i) {
            expect(D.entry(i, i)).toBe(0);
            for (let j = 0; j < i; ++j) {
                expect(D.entry(i, j)).toBe(D.entry(j, i));
            }
        }
    });

    it("KMeans agrees on both paths", () => {
        const X = random_data(60, 4);
        const { wasm, js } = both_paths(() => new KMeans(X, { K: 3, seed: 42 }).get_clusters());
        expect(wasm).toEqual(js);
    });

    // UMAP's SGD is chaotic — one ulp becomes O(10) within a few epochs, and how fast depends on the
    // graph, so comparing two full runs is meaningless no matter how few epochs are used. The kernel
    // is therefore tested directly against a JS transcription of `_optimize_layout` on a controlled
    // graph: same embedding, same edges, same pre-drawn negative samples, one epoch. Under those
    // conditions the two must agree to the last bit or two.
    //
    // The cases cover the branches that matter: edges skipped because they are not due this epoch,
    // and a varying number of negative samples per edge.
    for (const [label, vary] of [
        ["uniform schedule", false],
        ["mixed schedule with skipped edges", true],
    ]) {
        it(`UMAP SGD kernel matches the JS loop — ${label}`, () => {
            const R = new Randomizer(99);
            const n_points = 300;
            const n_edges = 900; // above WASM_MIN_UMAP_EDGES, so the kernel really runs
            expect(n_edges).toBeGreaterThanOrEqual(WASM_MIN_UMAP_EDGES);

            const dim = 2;
            const iter = 7;
            const a = 1.577;
            const b = 0.8951;
            const gamma = 1;
            const alpha = 0.9;

            const Y0 = Float64Array.from({ length: n_points * dim }, () => R.random * 4 - 2);
            const head = Int32Array.from({ length: n_edges }, () => R.random_int % n_points);
            const tail = Int32Array.from({ length: n_edges }, () => R.random_int % n_points);
            const eps = Float32Array.from({ length: n_edges }, () => (vary ? 1 + (R.random_int % 9) : 1));
            const eons = Float32Array.from({ length: n_edges }, () => (vary ? 1 + (R.random_int % 12) : 2));
            const epns = Float32Array.from({ length: n_edges }, () => (vary ? 1 + (R.random_int % 9) : 5));
            const eonns = Float32Array.from({ length: n_edges }, () => 1 + (R.random_int % 5));

            // Draw the negative samples exactly as UMAP does, in edge order.
            const negative = [];
            for (let i = 0; i < n_edges; ++i) {
                if (eons[i] > iter) continue;
                const n_neg = (iter - eonns[i]) / epns[i];
                for (let p = 0; p < n_neg; ++p) negative.push(R.random_int % n_points);
            }
            const neg = Int32Array.from(negative);

            const Y_wasm = Float64Array.from(Y0);
            const ok = wasmUmapOptimizeEpoch(
                Y_wasm,
                head,
                tail,
                eps,
                Float32Array.from(eons),
                epns,
                Float32Array.from(eonns),
                neg,
                dim,
                iter,
                a,
                b,
                gamma,
                alpha,
            );
            expect(ok, "kernel should have run").toBe(true);

            // JS transcription of UMAP._optimize_layout, replaying the same negative samples.
            const Y_js = Float64Array.from(Y0);
            const row = (i) => Y_js.subarray(i * dim, (i + 1) * dim);
            const sq = (p, q) => {
                let s = 0;
                for (let d = 0; d < dim; ++d) {
                    const t = p[d] - q[d];
                    s += t * t;
                }
                return s;
            };
            const clip = (x) => (x > 4 ? 4 : x < -4 ? -4 : x);
            let cursor = 0;
            for (let i = 0; i < n_edges; ++i) {
                if (eons[i] > iter) continue;
                const current = row(head[i]);
                const other = row(tail[i]);
                const dist = sq(current, other);
                if (dist > 0) {
                    const grad_coeff = (-2 * a * b * dist ** (b - 1)) / (a * dist ** b + 1);
                    for (let d = 0; d < dim; ++d) {
                        const g = clip(grad_coeff * (current[d] - other[d])) * alpha;
                        current[d] += g;
                        other[d] -= g;
                    }
                }
                const n_neg = (iter - eonns[i]) / epns[i];
                for (let p = 0; p < n_neg; ++p) {
                    const sample = row(neg[cursor++]);
                    const d2 = sq(current, sample);
                    if (d2 > 0) {
                        const grad_coeff = (2 * gamma * b) / ((0.01 + d2) * (a * d2 ** b + 1));
                        for (let d = 0; d < dim; ++d) {
                            const g = clip(grad_coeff * (current[d] - sample[d])) * alpha;
                            current[d] += g;
                            sample[d] -= g;
                        }
                    }
                }
            }

            expect(cursor, "kernel and reference must consume the same samples").toBe(neg.length);
            expect_close(Array.from(Y_wasm), Array.from(Y_js), `UMAP kernel (${label})`);
        });
    }

    // Iterative optimisers amplify the last-bit differences between the two accumulation orders, so
    // they are compared after a handful of steps rather than a full run.
    for (const [name, run] of [
        ["TSNE", (/** @type {number[][]} */ X) => new TSNE(X, { d: 2, seed: 42, perplexity: 10 }).transform(5)],
        ["SAMMON", (/** @type {number[][]} */ X) => new SAMMON(X, { d: 2, seed: 42 }).transform(5)],
        ["SMACOF", (/** @type {number[][]} */ X) => new SMACOF(X, { d: 2, seed: 42 }).transform(5)],
        ["SQDMDS", (/** @type {number[][]} */ X) => new SQDMDS(X, { d: 2, seed: 42 }).transform(5)],
        ["TriMap", (/** @type {number[][]} */ X) => new TriMap(X, { d: 2, seed: 42 }).transform(5)],
    ]) {
        it(`${name} agrees on both paths`, () => {
            const X = random_data(40, 5);
            const { wasm, js } = both_paths(() => run(X));
            expect_close(wasm, js, name);
        });
    }

    it("ISOMAP agrees on both paths", () => {
        // The JS side is a heap Dijkstra rather than a scalar mirror of the kernel, so this is a
        // stronger claim than the element-wise cases above: two different shortest-path
        // implementations have to produce the same geodesics.
        const X = random_data(40, 5);
        const { wasm, js } = both_paths(() => new ISOMAP(X, { d: 2, neighbors: 8, seed: 42 }).transform());
        expect_close(wasm, js, "ISOMAP");
    });

    it("ISOMAP declines the worker pool below the parallel threshold", () => {
        // `_dijkstra_parallel` is the pool entry point; under `WASM_MIN_PARALLEL_ROWS` it must bow
        // out so the single-threaded kernel handles the run.
        const isomap = new ISOMAP(random_data(20, 4), { d: 2, neighbors: 5, seed: 42 });
        const rows = 20;
        const out = new Float64Array(rows * rows);
        expect(isomap._dijkstra_parallel(new Int32Array(rows * 5), new Float64Array(rows * 5), out, rows, 5)).toBe(
            false,
        );
    });

    for (const kernel of ["gaussian", "flat"]) {
        it(`MeanShift agrees on both paths — ${kernel} kernel`, () => {
            const X = random_data(40, 4);
            const { wasm, js } = both_paths(() =>
                new MeanShift(X, { bandwidth: 0.4, kernel, max_iter: 5, seed: 42 }).get_cluster_list(),
            );
            expect(wasm).toEqual(js);
        });
    }

    it("MeanShift declines the worker pool below the parallel threshold", () => {
        const ms = new MeanShift(random_data(20, 4), { bandwidth: 0.4, max_iter: 3, seed: 42 });
        expect(ms._mean_shift_parallel(20, 4, true)).toBe(null);
    });

    it("MeanShift takes the metric-agnostic path for a non-euclidean metric", () => {
        // Every kernel here is euclidean-only, so a custom metric has to bypass WASM entirely and
        // still cluster — the branch that decides this is never taken by the default parameters.
        const X = random_data(40, 4);
        const { wasm, js } = both_paths(() =>
            new MeanShift(X, { bandwidth: 1.5, metric: manhattan, max_iter: 5, seed: 42 }).get_cluster_list(),
        );
        expect(wasm).toEqual(js);
    });
});
