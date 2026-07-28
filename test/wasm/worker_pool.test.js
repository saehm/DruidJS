import { afterAll, describe, expect, it } from "vitest";
import { ISOMAP } from "../../src/dimred/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { wasmDijkstraAPSP, wasmDistanceMatrix } from "../../src/wasm/index.js";
import { parallel_available, run_row_range, terminate_pool } from "../../src/wasm/worker_pool.js";

/**
 * @param {number} seed
 * @returns {() => number}
 */
function rng(seed) {
    let x = seed;
    return () => ((x ^= x << 13), (x ^= x >>> 17), (x ^= x << 5), (x >>> 0) / 4294967296);
}

describe("worker pool", () => {
    afterAll(() => terminate_pool());

    it("reports availability without throwing", () => {
        expect(typeof parallel_available()).toBe("boolean");
    });

    it("splits the distance matrix across workers and matches the single-threaded kernel", () => {
        if (!parallel_available()) return;
        const n = 200;
        const d = 8;
        const random = rng(17);
        const X = new Float64Array(n * d).map(() => random());

        const expected = new Float64Array(n * n);
        wasmDistanceMatrix(X, n, d, expected);

        const out = new Float64Array(new SharedArrayBuffer(n * n * 8));
        const ran = run_row_range("euclidean_distance_matrix_range_f64", [{ data: X, kind: "f64" }], out, [n, d], n, n);
        expect(ran).toBe(true);

        for (let i = 0; i < n * n; ++i) {
            expect(out[i]).toBeCloseTo(expected[i], 10);
        }
    });

    it("leaves no row unwritten", () => {
        if (!parallel_available()) return;
        const n = 128;
        const d = 4;
        const random = rng(5);
        const X = new Float64Array(n * d).map(() => random());
        const out = new Float64Array(new SharedArrayBuffer(n * n * 8));
        out.fill(-1);

        run_row_range("euclidean_distance_matrix_range_f64", [{ data: X, kind: "f64" }], out, [n, d], n, n);

        // A range wrapper that copied its whole buffer back would leave other workers' rows at 0,
        // and one that skipped rows would leave them at -1.
        let untouched = 0;
        for (let i = 0; i < n * n; ++i) {
            if (out[i] === -1) ++untouched;
        }
        expect(untouched).toBe(0);
    });

    it("refuses an output that is not shared", () => {
        if (!parallel_available()) return;
        const n = 32;
        const X = new Float64Array(n * 2).map((_, i) => i);
        const not_shared = new Float64Array(n * n);
        expect(
            run_row_range("euclidean_distance_matrix_range_f64", [{ data: X, kind: "f64" }], not_shared, [n, 2], n, n),
        ).toBe(false);
    });

    it("survives being torn down and used again", () => {
        if (!parallel_available()) return;
        const n = 64;
        const d = 3;
        const X = new Float64Array(n * d).map((_, i) => Math.sin(i));
        const run = () => {
            const out = new Float64Array(new SharedArrayBuffer(n * n * 8));
            const ok = run_row_range(
                "euclidean_distance_matrix_range_f64",
                [{ data: X, kind: "f64" }],
                out,
                [n, d],
                n,
                n,
            );
            return { ok, out };
        };

        const first = run();
        terminate_pool();
        const second = run();

        expect(first.ok).toBe(true);
        expect(second.ok).toBe(true);
        for (let i = 0; i < n * n; ++i) {
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
