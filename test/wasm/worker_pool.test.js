import { afterAll, describe, expect, it } from "vitest";
import { wasmMeanShiftStep } from "../../src/clustering/MeanShift.wasm.js";
import { wasmDijkstraAPSP } from "../../src/dimred/ISOMAP.wasm.js";
import { ISOMAP } from "../../src/dimred/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { parallel_available, run_row_range, terminate_pool } from "../../src/wasm/worker_pool.js";

/**
 * @param {number} seed
 * @returns {() => number}
 */
function rng(seed) {
    let x = seed;
    return () => ((x ^= x << 13), (x ^= x >>> 17), (x ^= x << 5), (x >>> 0) / 4294967296);
}

const BANDWIDTH = 0.7;
const D = 4;

/**
 * Splits one mean-shift step across the pool.
 *
 * The pool's own mechanics are what these tests are after, so they need some row-partitioned kernel
 * to drive it; mean shift is used because it is one of the two the library actually runs in
 * parallel, it writes every entry of its output (so an unwritten row is visible), and it reports a
 * per-worker maximum back through `returns`, which nothing else here covers.
 *
 * @param {Float64Array} points
 * @param {number} n
 * @param {Float64Array} [out] - Shared output; allocated if omitted.
 * @returns {{ ran: boolean; out: Float64Array; shifts: Float64Array }}
 */
function shift_step(points, n, out) {
    const shared = out ?? new Float64Array(new SharedArrayBuffer(n * D * 8));
    const shifts = new Float64Array(new SharedArrayBuffer(64));
    const inputs = /** @type {{ data: Float64Array | Int32Array; kind: "f64" | "i32" }[]} */ ([
        { data: points, kind: "f64" },
    ]);
    const ran = run_row_range("meanshift_step_range_f64", inputs, shared, [n, D, BANDWIDTH, 1], n, D, shifts);
    return { ran, out: shared, shifts };
}

describe("worker pool", () => {
    afterAll(() => terminate_pool());

    it("reports availability without throwing", () => {
        expect(typeof parallel_available()).toBe("boolean");
    });

    it("splits a step across workers and matches the single-threaded kernel", () => {
        if (!parallel_available()) return;
        const n = 200;
        const random = rng(17);
        const X = new Float64Array(n * D).map(() => random());

        const expected = new Float64Array(n * D);
        const expected_shift = wasmMeanShiftStep(X, expected, n, D, BANDWIDTH, true);

        const { ran, out, shifts } = shift_step(X, n);
        expect(ran).toBe(true);

        for (let i = 0; i < n * D; ++i) {
            expect(out[i]).toBeCloseTo(expected[i], 10);
        }
        // Each worker reports the largest shift over its own rows, so the whole-input maximum is
        // the maximum of the slots — a worker reporting over rows it did not own would exceed it.
        expect(Math.max(...shifts)).toBeCloseTo(/** @type {number} */ (expected_shift), 10);
    });

    it("leaves no row unwritten", () => {
        if (!parallel_available()) return;
        const n = 128;
        const random = rng(5);
        const X = new Float64Array(n * D).map(() => random());
        const out = new Float64Array(new SharedArrayBuffer(n * D * 8));
        out.fill(-1);

        shift_step(X, n, out);

        // A range wrapper that copied its whole buffer back would leave other workers' rows at 0,
        // and one that skipped rows would leave them at -1.
        let untouched = 0;
        for (let i = 0; i < n * D; ++i) {
            if (out[i] === -1) ++untouched;
        }
        expect(untouched).toBe(0);
    });

    it("refuses an output that is not shared", () => {
        if (!parallel_available()) return;
        const n = 32;
        const X = new Float64Array(n * D).map((_, i) => i);
        const not_shared = new Float64Array(n * D);
        expect(shift_step(X, n, not_shared).ran).toBe(false);
    });

    it("survives being torn down and used again", () => {
        if (!parallel_available()) return;
        const n = 64;
        const X = new Float64Array(n * D).map((_, i) => Math.sin(i));

        const first = shift_step(X, n);
        terminate_pool();
        const second = shift_step(X, n);

        expect(first.ran).toBe(true);
        expect(second.ran).toBe(true);
        for (let i = 0; i < n * D; ++i) {
            expect(second.out[i]).toBeCloseTo(first.out[i], 12);
        }
    });

    it("gives ISOMAP the same embedding as the single-threaded path", () => {
        const random = rng(99);
        const X = new Matrix(120, 6, () => random());
        // 120 rows is below the parallel threshold, so this also covers the fallback.
        const a = new ISOMAP(X, { d: 2, neighbors: 8, seed: 1212 }).transform();
        const b = new ISOMAP(X, { d: 2, neighbors: 8, seed: 1212 }).transform();
        expect(a.to2dArray()).toEqual(b.to2dArray());
    });
});

describe("parallel Dijkstra", () => {
    afterAll(() => terminate_pool());

    it("matches the single-threaded all-pairs result", () => {
        if (!parallel_available()) return;
        const n = 150;
        const k = 6;
        const random = rng(23);
        const neighbors = new Int32Array(n * k);
        const distances = new Float64Array(n * k);
        for (let i = 0; i < n; ++i) {
            for (let j = 0; j < k; ++j) {
                neighbors[i * k + j] = Math.floor(random() * n);
                distances[i * k + j] = 0.1 + random();
            }
        }

        const expected = new Float64Array(n * n).fill(Number.POSITIVE_INFINITY);
        wasmDijkstraAPSP(neighbors, distances, expected, n, k);

        const out = new Float64Array(new SharedArrayBuffer(n * n * 8));
        const ran = run_row_range(
            "dijkstra_apsp_range_f64",
            [
                { data: neighbors, kind: "i32" },
                { data: distances, kind: "f64" },
            ],
            out,
            [n, k],
            n,
            n,
        );
        expect(ran).toBe(true);
        for (let i = 0; i < n * n; ++i) {
            expect(out[i]).toBe(expected[i]);
        }
    });
});
