import { describe, expect, it } from "vitest";
import { LocalMAP, PaCMAP } from "../../src/dimred/index.js";
import { HNSW } from "../../src/knn/index.js";
import { euclidean } from "../../src/metrics/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("LocalMAP", () => {
    it("should reduce dimensionality to d=2 by default", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const result = new LocalMAP(data, { n_neighbors: 4, seed: 42 }).transform(50);

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 30000 }, () => {
        const data = generateTestData(20, 5);
        const result = new LocalMAP(data, { n_neighbors: 4, seed: 42 }).transform(50);
        expectValidValues(result, "LocalMAP");
    });

    it("generator should yield intermediate results", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const gen = new LocalMAP(data, { n_neighbors: 4, seed: 42 }).generator(10);

        let count = 0;
        for (const result of gen) {
            count++;
            expect(result).toHaveLength(15);
        }
        expect(count).toBe(10);
    });

    it("should work with static transform", { timeout: 30000 }, () => {
        const data = generateTestData(12, 5);
        const result = LocalMAP.transform(data, { n_neighbors: 4 });
        expect(result).toHaveLength(12);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with static generator", { timeout: 30000 }, () => {
        const data = generateTestData(12, 5);
        const result = LocalMAP.generator(data, { n_neighbors: 4 }).next().value;
        expect(result).toHaveLength(12);
    });

    it("should work with static transform_async", { timeout: 30000 }, async () => {
        const data = generateTestData(12, 5);
        const result = await LocalMAP.transform_async(data, { n_neighbors: 4 });
        expect(result).toHaveLength(12);
    });

    it("should respect d=3 parameter", { timeout: 30000 }, () => {
        const data = generateTestData(15, 6);
        const result = new LocalMAP(data, { n_neighbors: 4, d: 3, seed: 42 }).transform(30);
        expect(result[0]).toHaveLength(3);
    });

    it("should be reproducible for a given seed", { timeout: 30000 }, () => {
        const data = generateTestData(20, 4);
        const a = new LocalMAP(data, { n_neighbors: 4, seed: 7 }).transform(40);
        const b = new LocalMAP(data, { n_neighbors: 4, seed: 7 }).transform(40);
        expect(a).toEqual(b);
    });

    it("should accept a knn index", { timeout: 30000 }, () => {
        const data = generateTestData(30, 4);
        const knn = new HNSW(data, { metric: euclidean, seed: 42, ef: 100 });
        const result = new LocalMAP(data, { n_neighbors: 4, seed: 42, knn }).transform(30);
        expect(result).toHaveLength(30);
        expectValidValues(result, "LocalMAP");
    });

    it("should expose low_dist_thres as a parameter", () => {
        const data = generateTestData(15, 4);
        const dr = new LocalMAP(data, { n_neighbors: 4, low_dist_thres: 5 });
        expect(dr.parameter("low_dist_thres")).toBe(5);
        expect(new LocalMAP(data, { n_neighbors: 4 }).parameter("low_dist_thres")).toBe(10);
    });

    describe("fidelity to the reference implementation", () => {
        it("should be identical to PaCMAP through phases 1 and 2", { timeout: 30000 }, () => {
            // The reference switches gradients on `itr > phase_1 + phase_2`, so the first phase 3
            // step still runs plain PaCMAP: the two must agree for phase_1 + phase_2 + 1 steps.
            const data = generateTestData(30, 4);
            const opts = { n_neighbors: 4, seed: 42, num_iters: [10, 10, 20] };
            const steps = 10 + 10 + 1;

            const p = new PaCMAP(data, opts).transform(steps);
            const l = new LocalMAP(data, opts).transform(steps);
            expect(l).toEqual(p);
        });

        it("should diverge from PaCMAP once phase 3 is under way", { timeout: 30000 }, () => {
            const data = generateTestData(30, 4);
            const opts = { n_neighbors: 4, seed: 42, num_iters: [10, 10, 20] };

            const p = new PaCMAP(data, opts).transform(25);
            const l = new LocalMAP(data, opts).transform(25);
            expect(l).not.toEqual(p);
        });

        it("local nn gradient should match `pacmap_grad_nearby_recip_sqrt`", () => {
            const data = generateTestData(20, 4);
            const dr = new LocalMAP(data, { n_neighbors: 4, seed: 42, low_dist_thres: 8 });
            dr.init();

            const pairs = new Int32Array([0, 5, 1, 9, 7, 2, 11, 4, 15, 19]);
            const w_nn = 1.0;
            const nn_scale = 8 / 2;

            const [n, dim] = dr.Y.shape;
            const actual = new Float64Array(n * dim);
            dr._accumulate_gradients_local_nn(actual, pairs, w_nn, nn_scale);

            // Reference: w1 = w * (20/(10+d_ij)^2); w1 *= NN_coef_recip/sqrt(d_ij)
            const expected = new Float64Array(n * dim);
            for (let t = 0; t < pairs.length / 2; ++t) {
                const i = pairs[t * 2];
                const j = pairs[t * 2 + 1];
                const y_ij = new Float64Array(dim);
                let d_ij = 1.0;
                for (let k = 0; k < dim; ++k) {
                    y_ij[k] = dr.Y.entry(i, k) - dr.Y.entry(j, k);
                    d_ij += y_ij[k] ** 2;
                }
                let w1 = w_nn * (20 / (10 + d_ij) ** 2);
                w1 *= nn_scale / Math.sqrt(d_ij);
                for (let k = 0; k < dim; ++k) {
                    expected[i * dim + k] += w1 * y_ij[k];
                    expected[j * dim + k] -= w1 * y_ij[k];
                }
            }

            for (let i = 0; i < expected.length; ++i) {
                expect(actual[i]).toBeCloseTo(expected[i], 12);
            }
        });

        it("should only redraw further pairs within low_dist_thres", () => {
            const data = generateTestData(60, 4);
            const dr = new LocalMAP(data, { n_neighbors: 5, seed: 42, low_dist_thres: 0.05 });
            dr.init();

            const before = Int32Array.from(dr._fp_pairs);
            dr._resample_local_fp_pairs(0.05);

            const d = 2;
            const n_FP = dr._fp_pairs.length / 2 / 60;
            let changed = 0;
            for (let i = 0; i < 60; ++i) {
                for (let c = 0; c < n_FP; ++c) {
                    const t = i * n_FP + c;
                    const j = dr._fp_pairs[t * 2 + 1];
                    if (j === before[t * 2 + 1]) continue;
                    changed++;
                    // A replaced partner must be inside the threshold in the embedding.
                    let sq = 0;
                    for (let k = 0; k < d; ++k) {
                        const diff = dr.Y.values[i * d + k] - dr.Y.values[j * d + k];
                        sq += diff * diff;
                    }
                    expect(Math.sqrt(sq)).toBeLessThanOrEqual(0.05);
                }
            }
            expect(changed).toBeGreaterThan(0);
        });

        it("should keep further pairs distinct and off the neighbour list after redrawing", () => {
            const data = generateTestData(60, 4);
            const dr = new LocalMAP(data, { n_neighbors: 5, seed: 42, low_dist_thres: 100 });
            dr.init();
            dr._resample_local_fp_pairs(100);

            const n_FP = dr._fp_pairs.length / 2 / 60;
            for (let i = 0; i < 60; ++i) {
                const nn = new Set();
                for (let c = 0; c < 5; ++c) nn.add(dr._nn_pairs[(i * 5 + c) * 2 + 1]);
                const seen = new Set();
                for (let c = 0; c < n_FP; ++c) {
                    const j = dr._fp_pairs[(i * n_FP + c) * 2 + 1];
                    expect(j).not.toBe(i);
                    expect(nn.has(j)).toBe(false);
                    expect(seen.has(j)).toBe(false);
                    seen.add(j);
                }
            }
        });

        it("should redraw on exactly the reference's iterations", () => {
            // Reference: `if (itr > p1 + p2) and (itr % 10 == 0)`, where `itr` is the loop variable
            // -- the iteration that just ran, not the next one. Keying it off the post-increment
            // counter instead shifts every redraw one step early, which nothing else would catch.
            const data = generateTestData(40, 4);
            const dr = new LocalMAP(data, { n_neighbors: 4, seed: 42, num_iters: [5, 5, 30] });
            dr.init();

            const fired = [];
            const original = dr._resample_local_fp_pairs.bind(dr);
            dr._resample_local_fp_pairs = (thres) => {
                fired.push(dr._iter - 1); // the iteration that just completed
                return original(thres);
            };

            for (let i = 0; i < 40; ++i) dr.next();

            // p1 + p2 = 10, so: strictly after 10, every multiple of 10 -> 20 and 30.
            expect(fired).toEqual([20, 30]);
        });

        it("should keep the old partner when nothing falls within the threshold", () => {
            const data = generateTestData(40, 4);
            const dr = new LocalMAP(data, { n_neighbors: 4, seed: 42, low_dist_thres: 1e-12 });
            dr.init();

            const before = Int32Array.from(dr._fp_pairs);
            dr._resample_local_fp_pairs(1e-12);
            expect(Array.from(dr._fp_pairs)).toEqual(Array.from(before));
        });
    });
});
