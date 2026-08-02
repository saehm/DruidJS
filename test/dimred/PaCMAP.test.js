import { describe, expect, it } from "vitest";
import { PaCMAP } from "../../src/dimred/index.js";
import { HNSW, NaiveKNN, spatial_tree } from "../../src/knn/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { euclidean, manhattan } from "../../src/metrics/index.js";
import { setWasmEnabled } from "../../src/wasm/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("PaCMAP", () => {
    it("should reduce dimensionality to d=2 by default", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const pacmap = new PaCMAP(data, { n_neighbors: 4, seed: 42 });
        const result = pacmap.transform(50);

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should complete within timeout", { timeout: 30000 }, () => {
        const data = generateTestData(20, 5);
        const pacmap = new PaCMAP(data, { n_neighbors: 5, d: 2, seed: 42 });
        const result = pacmap.transform(50);

        expect(result).toHaveLength(20);
        expect(result[0]).toHaveLength(2);
    });

    it("generator should yield intermediate results", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const pacmap = new PaCMAP(data, { n_neighbors: 5, d: 2, seed: 42 });
        const gen = pacmap.generator(10);

        let count = 0;
        for (const result of gen) {
            count++;
            expect(result).toHaveLength(15);
            expect(result[0]).toHaveLength(2);
        }
        expect(count).toBe(10);
    });

    it("should produce valid values", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const pacmap = new PaCMAP(data, { n_neighbors: 4, d: 2, seed: 42 });
        const result = pacmap.transform(50);
        expectValidValues(result, "PaCMAP");
    });

    it("should work with static transform", { timeout: 30000 }, () => {
        const data = generateTestData(12, 5);
        const result = PaCMAP.transform(data, { n_neighbors: 4 });
        expect(result).toHaveLength(12);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with static generator", { timeout: 30000 }, () => {
        const data = generateTestData(12, 5);
        const gen = PaCMAP.generator(data, { n_neighbors: 4 });
        const result = gen.next().value;
        expect(result).toHaveLength(12);
    });

    it("should work with static transform_async", { timeout: 30000 }, async () => {
        const data = generateTestData(12, 5);
        const result = await PaCMAP.transform_async(data, { n_neighbors: 4 });
        expect(result).toHaveLength(12);
    });

    it("should respect d=3 parameter", { timeout: 30000 }, () => {
        const data = generateTestData(15, 6);
        const pacmap = new PaCMAP(data, { n_neighbors: 4, d: 3, seed: 42 });
        const result = pacmap.transform(30);

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(3);
    });

    it("should throw if n_neighbors >= N", () => {
        const data = generateTestData(10, 4);
        expect(() => new PaCMAP(data, { n_neighbors: 10 })).toThrow();
    });

    it("should be reproducible for a given seed", { timeout: 30000 }, () => {
        const data = generateTestData(20, 4);
        const a = new PaCMAP(data, { n_neighbors: 4, seed: 7 }).transform(40);
        const b = new PaCMAP(data, { n_neighbors: 4, seed: 7 }).transform(40);
        expect(a).toEqual(b);
    });

    it("should accept an exact knn index and match the default", { timeout: 30000 }, () => {
        const data = generateTestData(30, 4);
        const expected = new PaCMAP(data, { n_neighbors: 4, seed: 42 }).transform(30);

        for (const knn of [
            spatial_tree(data, { metric: euclidean, seed: 42 }),
            new NaiveKNN(data, { metric: euclidean, seed: 42 }),
        ]) {
            const actual = new PaCMAP(data, { n_neighbors: 4, seed: 42, knn }).transform(30);
            expect(actual).toEqual(expected);
        }
    });

    it("should accept an approximate knn index", { timeout: 30000 }, () => {
        const data = generateTestData(30, 4);
        const knn = new HNSW(data, { metric: euclidean, seed: 42, ef: 100 });
        const result = new PaCMAP(data, { n_neighbors: 4, seed: 42, knn }).transform(30);
        expect(result).toHaveLength(30);
        expectValidValues(result, "PaCMAP");
    });

    it("should agree between the WASM and JS neighbour search", { timeout: 30000 }, () => {
        const data = generateTestData(30, 4);
        const expected = new PaCMAP(data, { n_neighbors: 4, seed: 42 }).transform(30);
        setWasmEnabled(false);
        try {
            const actual = new PaCMAP(data, { n_neighbors: 4, seed: 42 }).transform(30);
            expect(actual).toEqual(expected);
        } finally {
            setWasmEnabled(true);
        }
    });

    it("should support a non-euclidean metric", { timeout: 30000 }, () => {
        const data = generateTestData(20, 4);
        const result = new PaCMAP(data, { n_neighbors: 4, seed: 42, metric: manhattan }).transform(30);
        expect(result).toHaveLength(20);
        expectValidValues(result, "PaCMAP");
    });

    describe("fidelity to the reference implementation", () => {
        /** Reference `pacmap_grad`, transcribed from YingfanWang/PaCMAP. */
        function reference_grad(Y, pairs, kind, w) {
            const [n, dim] = Y.shape;
            const grad = new Float64Array(n * dim);
            for (let t = 0; t < pairs.length / 2; ++t) {
                const i = pairs[t * 2];
                const j = pairs[t * 2 + 1];
                const y_ij = new Float64Array(dim);
                let d_ij = 1.0;
                for (let k = 0; k < dim; ++k) {
                    y_ij[k] = Y.entry(i, k) - Y.entry(j, k);
                    d_ij += y_ij[k] ** 2;
                }
                if (kind === "nn" || kind === "mn") {
                    const c = kind === "nn" ? 10 : 10000;
                    const w1 = w * ((2 * c) / (c + d_ij) ** 2);
                    for (let k = 0; k < dim; ++k) {
                        grad[i * dim + k] += w1 * y_ij[k];
                        grad[j * dim + k] -= w1 * y_ij[k];
                    }
                } else {
                    const w1 = (w * 2) / (1 + d_ij) ** 2;
                    for (let k = 0; k < dim; ++k) {
                        grad[i * dim + k] -= w1 * y_ij[k];
                        grad[j * dim + k] += w1 * y_ij[k];
                    }
                }
            }
            return grad;
        }

        const cases = [
            ["nn", 10, false, 2.0],
            ["mn", 10000, false, 3.0],
            ["fp", 1, true, 1.0],
        ];

        for (const [kind, a, repulsive, w] of cases) {
            it(`${kind} gradient should match the reference formula`, () => {
                const data = generateTestData(20, 4);
                const dr = new PaCMAP(data, { n_neighbors: 4, seed: 42 });
                dr.init();
                const pairs = new Int32Array([0, 5, 1, 9, 3, 3, 7, 2, 11, 4, 15, 19]);

                const actual = new Float64Array(dr.Y.shape[0] * dr.Y.shape[1]);
                dr._accumulate_gradients(actual, pairs, w, a, repulsive);
                const expected = reference_grad(dr.Y, pairs, kind, w);

                for (let i = 0; i < expected.length; ++i) {
                    expect(actual[i]).toBeCloseTo(expected[i], 12);
                }
            });
        }

        it("should scale neighbour distances by the local density", () => {
            const data = generateTestData(40, 4);
            const dr = new PaCMAP(data, { n_neighbors: 5, seed: 42 });
            dr.init();

            const n_candidates = Math.min(5 + 50, 40 - 1);
            const { dist, idx } = dr._find_candidates(dr._X_knn, n_candidates);

            // sigma is the mean over the 4th to 6th neighbour, floored, as in `generate_pair`.
            const sig = [];
            for (let i = 0; i < 40; ++i) {
                let sum = 0;
                for (let c = 3; c < 6; ++c) sum += dist[i * n_candidates + c];
                sig.push(Math.max(sum / 3, 1e-10));
            }

            for (let i = 0; i < 40; ++i) {
                const scaled = [];
                for (let c = 0; c < n_candidates; ++c) {
                    const j = idx[i * n_candidates + c];
                    const dv = dist[i * n_candidates + c];
                    scaled.push([(dv * dv) / (sig[i] * sig[j]), j]);
                }
                scaled.sort((p, q) => p[0] - q[0]);
                for (let c = 0; c < 5; ++c) {
                    expect(dr._nn_pairs[(i * 5 + c) * 2]).toBe(i);
                    expect(dr._nn_pairs[(i * 5 + c) * 2 + 1]).toBe(scaled[c][1]);
                }
            }
        });

        it("should select neighbours that differ from an unscaled knn graph", () => {
            // If the two ever agree on every point the rescaling has silently stopped working.
            const data = generateTestData(60, 5);
            const dr = new PaCMAP(data, { n_neighbors: 5, seed: 42 });
            dr.init();
            const n_candidates = Math.min(5 + 50, 60 - 1);
            const { idx } = dr._find_candidates(dr._X_knn, n_candidates);

            let differs = 0;
            for (let i = 0; i < 60; ++i) {
                for (let c = 0; c < 5; ++c) {
                    if (dr._nn_pairs[(i * 5 + c) * 2 + 1] !== idx[i * n_candidates + c]) differs++;
                }
            }
            expect(differs).toBeGreaterThan(0);
        });

        it("should normalize the input as `preprocess_X` does", () => {
            const data = generateTestData(30, 4).map((row) => row.map((v) => v * 1000 + 500));
            const dr = new PaCMAP(data, { n_neighbors: 4, seed: 42 });
            dr.init();

            // Reference: X -= min(X); X /= max(X); X -= mean(X, axis=0). The per-column centring
            // comes last, so the global range is only unit *before* it.
            const flat = data.flat();
            const xmin = Math.min(...flat);
            const xmax = Math.max(...flat.map((v) => v - xmin));
            const expected = data.map((row) => row.map((v) => (v - xmin) / xmax));
            for (let j = 0; j < 4; ++j) {
                let mean = 0;
                for (let i = 0; i < 30; ++i) mean += expected[i][j];
                mean /= 30;
                for (let i = 0; i < 30; ++i) expected[i][j] -= mean;
            }

            for (let i = 0; i < 30; ++i) {
                for (let j = 0; j < 4; ++j) {
                    expect(dr._X_knn.entry(i, j)).toBeCloseTo(expected[i][j], 12);
                }
            }
        });

        it("should scale the input to the loss constants, not the raw units", () => {
            // Same shape, different units: normalization makes the layout invariant to it.
            const base = generateTestData(25, 4);
            const scaled = base.map((row) => row.map((v) => v * 250));
            const a = new PaCMAP(base, { n_neighbors: 4, seed: 42 }).transform(40);
            const b = new PaCMAP(scaled, { n_neighbors: 4, seed: 42 }).transform(40);
            for (let i = 0; i < a.length; ++i) {
                for (let k = 0; k < 2; ++k) expect(a[i][k]).toBeCloseTo(b[i][k], 8);
            }
        });

        it("should draw mid-near pairs from outside the neighbour list on average", () => {
            const data = generateTestData(80, 5);
            const dr = new PaCMAP(data, { n_neighbors: 5, seed: 42 });
            dr.init();

            // Reference MN pairs are the 2nd closest of six *global* draws, so they must be
            // typically farther than the NN pairs but nearer than a uniform random point.
            const mean_of = (pairs) => {
                let sum = 0;
                for (let t = 0; t < pairs.length / 2; ++t) {
                    sum += euclidean(dr._X_knn.row(pairs[t * 2]), dr._X_knn.row(pairs[t * 2 + 1]));
                }
                return sum / (pairs.length / 2);
            };
            const nn = mean_of(dr._nn_pairs);
            const mn = mean_of(dr._mn_pairs);
            const fp = mean_of(dr._fp_pairs);
            expect(mn).toBeGreaterThan(nn);
            expect(mn).toBeLessThan(fp);
        });

        it("should never pair a point with itself or draw a further pair from its neighbours", () => {
            const data = generateTestData(40, 4);
            const dr = new PaCMAP(data, { n_neighbors: 5, seed: 42 });
            dr.init();

            for (const pairs of [dr._nn_pairs, dr._mn_pairs, dr._fp_pairs]) {
                for (let t = 0; t < pairs.length / 2; ++t) {
                    expect(pairs[t * 2]).not.toBe(pairs[t * 2 + 1]);
                }
            }
            for (let i = 0; i < 40; ++i) {
                const nn = new Set();
                for (let c = 0; c < 5; ++c) nn.add(dr._nn_pairs[(i * 5 + c) * 2 + 1]);
                const n_FP = dr._fp_pairs.length / 2 / 40;
                const seen = new Set();
                for (let c = 0; c < n_FP; ++c) {
                    const j = dr._fp_pairs[(i * n_FP + c) * 2 + 1];
                    expect(nn.has(j)).toBe(false);
                    expect(seen.has(j)).toBe(false);
                    seen.add(j);
                }
            }
        });

        it("should follow the reference weight schedule", () => {
            const data = generateTestData(15, 4);
            const dr = new PaCMAP(data, { n_neighbors: 4, seed: 42 });
            const [p1, p2] = dr.parameter("num_iters");

            expect(dr._get_weights(0)).toEqual({ w_nn: 2.0, w_mn: 1000.0, w_fp: 1.0 });
            expect(dr._get_weights(p1 - 1).w_mn).toBeCloseTo(1000 * (1 / p1) + 3 * (1 - 1 / p1), 10);
            expect(dr._get_weights(p1)).toEqual({ w_nn: 3.0, w_mn: 3.0, w_fp: 1.0 });
            expect(dr._get_weights(p1 + p2 - 1)).toEqual({ w_nn: 3.0, w_mn: 3.0, w_fp: 1.0 });
            expect(dr._get_weights(p1 + p2)).toEqual({ w_nn: 1.0, w_mn: 0.0, w_fp: 1.0 });
        });

        it("should apply Adam with the reference bias correction", () => {
            const data = generateTestData(15, 4);
            const dr = new PaCMAP(data, { n_neighbors: 4, seed: 42, lr: 0.5 });
            dr.init();

            const Y0 = dr.Y.clone();
            const total = dr.Y.shape[0] * dr.Y.shape[1];
            const grad = new Float64Array(total).map((_, i) => Math.sin(i) * 0.3);
            dr._adam_update(grad);

            // First step: m = 0.1*g, v = 0.001*g², lr_t = lr*sqrt(1-0.999)/(1-0.9)
            const lr_t = (0.5 * Math.sqrt(1 - 0.999)) / (1 - 0.9);
            for (let i = 0; i < total; ++i) {
                const m = 0.1 * grad[i];
                const v = 0.001 * grad[i] ** 2;
                const expected = Y0.values[i] - (lr_t * m) / (Math.sqrt(v) + 1e-7);
                expect(dr.Y.values[i]).toBeCloseTo(expected, 12);
            }
        });

        it("should initialize from PCA scaled by 0.01", () => {
            const data = generateTestData(30, 6);
            const dr = new PaCMAP(data, { n_neighbors: 4, seed: 42 });
            dr.init();

            // Y is 0.01 * the PCA of the preprocessed input, so its spread is small and centred.
            const [N, d] = dr.Y.shape;
            for (let k = 0; k < d; ++k) {
                let mean = 0;
                for (let i = 0; i < N; ++i) mean += dr.Y.entry(i, k);
                expect(mean / N).toBeCloseTo(0, 8);
            }
            expect(Math.max(...dr.Y.values.map(Math.abs))).toBeLessThan(0.1);
        });

        it("should reduce inputs wider than 100 dimensions before searching", { timeout: 30000 }, () => {
            const data = generateTestData(20, 120);
            const dr = new PaCMAP(data, { n_neighbors: 4, seed: 42 });
            dr.init();
            expect(dr._X_knn.shape[1]).toBe(100);

            const off = new PaCMAP(data, { n_neighbors: 4, seed: 42, apply_pca: false });
            off.init();
            expect(off._X_knn.shape[1]).toBe(120);
        });
    });

    it("should accept a Matrix and return a Matrix", { timeout: 30000 }, () => {
        const data = Matrix.from(generateTestData(20, 4));
        const result = new PaCMAP(data, { n_neighbors: 4, seed: 42 }).transform(20);
        expect(result).toBeInstanceOf(Matrix);
        expect(result.shape).toEqual([20, 2]);
    });
});
