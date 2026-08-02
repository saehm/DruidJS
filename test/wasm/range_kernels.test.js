import { describe, expect, it } from "vitest";
import { wasmMeanShiftStep, wasmMeanShiftStepRange } from "../../src/clustering/MeanShift.wasm.js";
import { wasmDijkstraAPSP, wasmDijkstraAPSPRange } from "../../src/dimred/ISOMAP.wasm.js";
import { wasmSqdmdsFillGrads, wasmSqdmdsFillGradsRange } from "../../src/dimred/SQDMDS.wasm.js";
import { wasmDistanceMatrix, wasmKnnBlock } from "../../src/matrix/distance_matrix.wasm.js";
import { wasmMatMul, wasmMatMulRange } from "../../src/matrix/Matrix.wasm.js";
import { Randomizer } from "../../src/util/randomizer.js";
import { isWasmAvailable } from "../../src/wasm/index.js";

/**
 * The `*Range` kernels are the halves the worker pool calls: each worker runs one row range against
 * a shared output buffer. Of the five, only `meanshift` and `dijkstra` are run in parallel by
 * anything today, so for the rest this file is the only thing exercising the row-offset arithmetic.
 *
 * Two properties are asserted for the kernels whose ranges partition the output rows — `matmul`,
 * `knn_block`, `meanshift` and `dijkstra`:
 *
 * 1. **Partition equivalence.** Splitting the rows into disjoint ranges and running each produces
 *    exactly what one whole-range call produces. An off-by-one in a range bound shows up here.
 * 2. **No clobbering.** A single range writes its own rows and *only* its own rows. Each wrapper
 *    carries a comment promising this ("so concurrent workers sharing `C_val` do not clobber each
 *    other's rows") because the kernel computes into a private buffer whose untouched rows are
 *    zero; copying the whole buffer back would erase every other worker's output. That is a silent
 *    data race, invisible single-threaded, so it is pinned here directly.
 *
 * The equality is exact, not approximate: both paths run the same kernel over the same lanes, so
 * the only difference is which rows get copied out.
 *
 * `sqdmds_fill_grads` is the exception. Its ranges partition *quartets*, and a quartet scatters
 * gradient into the four rows it names, so ranges overlap in the output and no slice contract
 * applies. It is checked for additivity instead — see that block.
 */

const SENTINEL = -12345.5;

/** @param {number} n @param {number} d @param {number} [seed] */
function random_flat(n, d, seed = 1212) {
    const R = new Randomizer(seed);
    return Float64Array.from({ length: n * d }, () => R.random);
}

/**
 * Splits `[0, rows)` into contiguous ranges, including an empty one and a single-row one so the
 * degenerate bounds are covered too.
 *
 * @param {number} rows
 * @returns {[number, number][]}
 */
function ranges(rows) {
    const mid = rows >> 1;
    return [
        [0, 1],
        [1, mid],
        [mid, mid], // empty range: must write nothing at all
        [mid, rows],
    ];
}

describe.skipIf(!isWasmAvailable())("range kernels", () => {
    describe("wasmMatMulRange", () => {
        const rows_A = 12;
        const cols_A = 7;
        const cols_B = 5;
        const A = random_flat(rows_A, cols_A, 1);
        const B = random_flat(cols_A, cols_B, 2);

        const full = new Float64Array(rows_A * cols_B);
        wasmMatMul(A, rows_A, cols_A, B, cols_B, full);

        it("should reproduce the whole product from disjoint row ranges", () => {
            const split = new Float64Array(rows_A * cols_B).fill(SENTINEL);
            for (const [start, end] of ranges(rows_A)) {
                expect(wasmMatMulRange(A, rows_A, cols_A, B, cols_B, split, start, end)).toBe(true);
            }
            expect(Array.from(split)).toEqual(Array.from(full));
        });

        it("should leave rows outside its range untouched", () => {
            const out = new Float64Array(rows_A * cols_B).fill(SENTINEL);
            const start = 3;
            const end = 8;
            wasmMatMulRange(A, rows_A, cols_A, B, cols_B, out, start, end);

            for (let i = 0; i < rows_A; ++i) {
                for (let j = 0; j < cols_B; ++j) {
                    const at = i * cols_B + j;
                    const expected = i >= start && i < end ? full[at] : SENTINEL;
                    expect(out[at], `row ${i} col ${j}`).toBe(expected);
                }
            }
        });
    });

    describe("wasmKnnBlock", () => {
        const n = 14;
        const d = 4;
        const k = 5;
        const row_len = 2 * k;
        const X = random_flat(n, d, 3);

        const full = new Float64Array(n * row_len);
        wasmKnnBlock(X, n, d, k, full, 0, n);

        it("should reproduce the whole result from disjoint row ranges", () => {
            const split = new Float64Array(n * row_len).fill(SENTINEL);
            for (const [start, end] of ranges(n)) {
                expect(wasmKnnBlock(X, n, d, k, split, start, end)).toBe(true);
            }
            expect(Array.from(split)).toEqual(Array.from(full));
        });

        it("should leave rows outside its range untouched", () => {
            const out = new Float64Array(n * row_len).fill(SENTINEL);
            const start = 5;
            const end = 9;
            wasmKnnBlock(X, n, d, k, out, start, end);

            for (let i = 0; i < n; ++i) {
                for (let c = 0; c < row_len; ++c) {
                    const at = i * row_len + c;
                    const expected = i >= start && i < end ? full[at] : SENTINEL;
                    expect(out[at], `row ${i} slot ${c}`).toBe(expected);
                }
            }
        });

        it("should return the exact k nearest, ascending, excluding the row itself", () => {
            // Selection happens inside the kernel, so this is the only place it is checked against
            // a full sort of every distance.
            const D = new Float64Array(n * n);
            wasmDistanceMatrix(X, n, d, D);

            for (let i = 0; i < n; ++i) {
                const exact = [];
                for (let j = 0; j < n; ++j) if (j !== i) exact.push([D[i * n + j], j]);
                exact.sort((a, b) => a[0] - b[0]);

                for (let c = 0; c < k; ++c) {
                    const got_dist = full[i * row_len + c];
                    const got_idx = full[i * row_len + k + c];
                    expect(got_idx, `row ${i} slot ${c}`).not.toBe(i);
                    expect(got_dist, `row ${i} slot ${c}`).toBeCloseTo(exact[c][0], 12);
                    // Ties may order differently, so verify the index really is at that distance.
                    expect(D[i * n + got_idx]).toBeCloseTo(exact[c][0], 12);
                    if (c > 0) expect(got_dist).toBeGreaterThanOrEqual(full[i * row_len + c - 1]);
                }
            }
        });

        it("should pad with -1 when fewer than k other points exist", () => {
            const small_n = 4;
            const big_k = 6;
            const Xs = random_flat(small_n, d, 9);
            const out = new Float64Array(small_n * 2 * big_k);
            wasmKnnBlock(Xs, small_n, d, big_k, out, 0, small_n);

            for (let i = 0; i < small_n; ++i) {
                for (let c = small_n - 1; c < big_k; ++c) {
                    expect(out[i * 2 * big_k + big_k + c], `row ${i} slot ${c}`).toBe(-1);
                }
            }
        });
    });

    describe("wasmMeanShiftStepRange", () => {
        const n = 16;
        const d = 3;
        const bandwidth = 0.5;
        const points = random_flat(n, d, 4);

        for (const use_gaussian of [true, false]) {
            describe(use_gaussian ? "gaussian kernel" : "flat kernel", () => {
                const full = new Float64Array(n * d);
                const full_shift = wasmMeanShiftStep(points, full, n, d, bandwidth, use_gaussian);

                it("should reproduce the whole step from disjoint point ranges", () => {
                    const split = new Float64Array(n * d).fill(SENTINEL);
                    for (const [start, end] of ranges(n)) {
                        expect(
                            wasmMeanShiftStepRange(points, split, n, d, bandwidth, use_gaussian, start, end),
                        ).not.toBe(null);
                    }
                    expect(Array.from(split)).toEqual(Array.from(full));
                });

                it("should report the max shift over its own range only", () => {
                    // Every range's max shift is bounded by the global one, and the range covering
                    // everything must reproduce it exactly.
                    const scratch = new Float64Array(n * d);
                    const whole = wasmMeanShiftStepRange(points, scratch, n, d, bandwidth, use_gaussian, 0, n);
                    expect(whole).toBe(full_shift);

                    for (const [start, end] of ranges(n)) {
                        if (start === end) continue;
                        const shift = wasmMeanShiftStepRange(
                            points,
                            new Float64Array(n * d),
                            n,
                            d,
                            bandwidth,
                            use_gaussian,
                            start,
                            end,
                        );
                        expect(shift).toBeLessThanOrEqual(/** @type {number} */ (full_shift));
                    }
                });

                it("should leave points outside its range untouched", () => {
                    const out = new Float64Array(n * d).fill(SENTINEL);
                    const start = 4;
                    const end = 10;
                    wasmMeanShiftStepRange(points, out, n, d, bandwidth, use_gaussian, start, end);

                    for (let i = 0; i < n; ++i) {
                        for (let j = 0; j < d; ++j) {
                            const at = i * d + j;
                            const expected = i >= start && i < end ? full[at] : SENTINEL;
                            expect(out[at], `point ${i} dim ${j}`).toBe(expected);
                        }
                    }
                });
            });
        }
    });

    describe("wasmDijkstraAPSPRange", () => {
        // A ring graph: every node is joined to its two neighbours, so the geodesic between any two
        // nodes is the shorter way round and no node is unreachable.
        const n = 12;
        const k = 2;
        const neighbors = new Int32Array(n * k);
        const distances = new Float64Array(n * k);
        for (let i = 0; i < n; ++i) {
            neighbors[i * k] = (i + 1) % n;
            neighbors[i * k + 1] = (i - 1 + n) % n;
            distances[i * k] = 1;
            distances[i * k + 1] = 1;
        }

        const full = new Float64Array(n * n);
        wasmDijkstraAPSP(neighbors, distances, full, n, k);

        it("should compute the ring's geodesics", () => {
            // Sanity check on the fixture itself, so a broken kernel cannot be masked by both paths
            // agreeing on nonsense.
            for (let i = 0; i < n; ++i) {
                for (let j = 0; j < n; ++j) {
                    const hops = Math.abs(i - j);
                    expect(full[i * n + j], `${i} -> ${j}`).toBeCloseTo(Math.min(hops, n - hops), 10);
                }
            }
        });

        it("should reproduce the whole APSP matrix from disjoint source ranges", () => {
            const split = new Float64Array(n * n).fill(SENTINEL);
            for (const [start, end] of ranges(n)) {
                expect(wasmDijkstraAPSPRange(neighbors, distances, split, n, k, start, end)).toBe(true);
            }
            expect(Array.from(split)).toEqual(Array.from(full));
        });

        it("should leave source rows outside its range untouched", () => {
            const out = new Float64Array(n * n).fill(SENTINEL);
            const start = 2;
            const end = 7;
            wasmDijkstraAPSPRange(neighbors, distances, out, n, k, start, end);

            for (let i = 0; i < n; ++i) {
                for (let j = 0; j < n; ++j) {
                    const at = i * n + j;
                    const expected = i >= start && i < end ? full[at] : SENTINEL;
                    expect(out[at], `source ${i} target ${j}`).toBe(expected);
                }
            }
        });
    });

    describe("wasmSqdmdsFillGradsRange", () => {
        // The odd one out: a quartet range scatters gradient into the four rows it names, so ranges
        // overlap in the output and cannot share a buffer. Each call returns the whole gradient for
        // its own quartets, and partials are summed — so the properties to pin are additivity and a
        // zero result for an empty range, not the disjoint-slice contract the others hold to.
        const n = 16;
        const d_ld = 2;
        const d_hd = 4;
        const num_q = 4;

        const rnd = new Randomizer(1212);
        const Y = Float64Array.from({ length: n * d_ld }, () => rnd.random);
        const X = Float64Array.from({ length: n * d_hd }, () => rnd.random);
        const quartets = new Uint32Array(num_q * 4);
        for (let q = 0; q < num_q; ++q) {
            const picked = [];
            while (picked.length < 4) {
                const c = Math.floor(rnd.random * n);
                if (!picked.includes(c)) picked.push(c);
            }
            quartets.set(picked, q * 4);
        }

        /** Runs a range and returns the gradient it produced. */
        function grads_for(start_q, end_q) {
            const g = new Float64Array(n * d_ld);
            expect(wasmSqdmdsFillGradsRange(Y, X, quartets, g, n, d_ld, d_hd, false, false, start_q, end_q)).toBe(true);
            return g;
        }

        const full = new Float64Array(n * d_ld);
        wasmSqdmdsFillGrads(Y, X, quartets, full, n, d_ld, d_hd, false, false);

        it("should match the non-range kernel when covering every quartet", () => {
            expect(Array.from(grads_for(0, num_q))).toEqual(Array.from(full));
        });

        it("should produce nothing for an empty range", () => {
            // The kernel accumulates without zeroing, so the wrapper has to clear its buffer. A
            // recycled allocation would otherwise surface here as leftovers from an earlier call.
            expect(Array.from(grads_for(2, 2))).toEqual(new Array(n * d_ld).fill(0));
        });

        it("should sum partial gradients to the whole gradient", () => {
            const a = grads_for(0, 2);
            const b = grads_for(2, num_q);
            for (let i = 0; i < full.length; ++i) {
                expect(a[i] + b[i], `entry ${i}`).toBeCloseTo(full[i], 12);
            }
        });

        it("should not depend on what ran before it", () => {
            // Same range, once from a clean slate and once after a different range has been through
            // the allocator: identical, or the buffer is leaking state between calls.
            const first = grads_for(0, 2);
            grads_for(2, num_q);
            expect(Array.from(grads_for(0, 2))).toEqual(Array.from(first));
        });
    });
});
