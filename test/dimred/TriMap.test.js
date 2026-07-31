import { describe, expect, it } from "vitest";
import { TriMap } from "../../src/dimred/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("TriMap", () => {
    it("should complete within timeout", { timeout: 30000 }, () => {
        const data = generateTestData(20, 5);
        const trimap = new TriMap(data, { d: 2, c: 2, seed: 42 });
        const result = trimap.transform(50);

        expect(result).toHaveLength(20);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const trimap = new TriMap(data, { d: 2, c: 2, seed: 42 });
        const result = trimap.transform(30);
        expectValidValues(result, "TriMap");
    });
    it("should work with static transform", () => {
        const data = generateTestData(60, 4);
        const result = TriMap.transform(data, { d: 2 });
        expect(result).toHaveLength(60);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with static generator", () => {
        const data = generateTestData(60, 4);
        const gen = TriMap.generator(data, { d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(60);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(60, 4);
        const result = await TriMap.transform_async(data, { d: 2 });
        expect(result).toHaveLength(60);
    });

    // Cross-checked against trimap 1.1.5: sig, the log-similarity matrix P and the triplet weights
    // agree with the reference to machine precision on identical neighbour data.
    describe("trimap 1.1.5 semantics", () => {
        const data = generateTestData(80, 5);
        const dr = new TriMap(data, { d: 2, n_inliers: 6, n_outliers: 3, n_random: 2, seed: 42 });
        dr.init();

        it("should keep similarities in log space", () => {
            // find_p returns -d^2 / (sig_i * sig_j) rather than its exponential, so entries are
            // <= 0 and large distances degrade linearly instead of underflowing to zero.
            const knn_distances = Matrix.from([
                [1, 2],
                [3, 400],
            ]);
            const nbrs = Matrix.from([
                [1, 0],
                [0, 1],
            ]);
            const sig = new Float64Array([1, 1]);
            // @ts-ignore - private
            const P = dr._find_p(knn_distances, sig, nbrs);

            expect(P.entry(0, 0)).toBeCloseTo(-1, 12);
            expect(P.entry(0, 1)).toBeCloseTo(-4, 12);
            expect(P.entry(1, 0)).toBeCloseTo(-9, 12);
            // exp(-160000) would have flushed to 0 and lost the ordering entirely.
            expect(P.entry(1, 1)).toBeCloseTo(-160000, 6);
        });

        it("should scale the PCA initialisation down", () => {
            // A raw PCA projection would dwarf the `1 +` offset in the loss and flatten gradients.
            const raw = generateTestData(80, 5);
            const scaled = new TriMap(raw, { d: 2, seed: 42 });
            scaled.init();
            const spread = Math.max(...scaled.Y.values.map(Math.abs));
            const dataSpread = Math.max(...raw.flat().map(Math.abs));
            expect(spread).toBeLessThan(dataSpread);
        });

        it("should produce non-negative weights whose minimum is zero", () => {
            // Weights are shifted by their minimum and passed through tempered_log(1 + w, 0.5),
            // so the smallest weight maps exactly to 0 and none can be negative.
            const weights = Array.from(dr.weights);
            expect(Math.min(...weights)).toBeCloseTo(0, 12);
            expect(weights.every((w) => w >= 0)).toBe(true);
            expect(weights.every((w) => Number.isFinite(w))).toBe(true);
        });

        it("should hold the learning rate fixed under delta-bar-delta", () => {
            const lr0 = dr.lr;
            // @ts-ignore - private
            dr._next(0);
            // @ts-ignore - private
            dr._next(1);
            expect(dr.lr).toBe(lr0);
        });

        it("should separate well-separated clusters", { timeout: 30000 }, () => {
            const centers = [
                [0, 0, 0, 0],
                [12, 12, 0, 0],
                [0, 0, 12, 12],
            ];
            const X = [];
            const labels = [];
            for (let c = 0; c < 3; ++c) {
                for (let i = 0; i < 40; ++i) {
                    X.push(centers[c].map((v, k) => v + Math.sin(i * 3.1 + k * 1.7 + c)));
                    labels.push(c);
                }
            }
            const Y = new TriMap(X, { d: 2, seed: 1212 }).transform(200);

            let agree = 0;
            for (let i = 0; i < X.length; ++i) {
                let best = -1;
                let bestD = Infinity;
                for (let j = 0; j < X.length; ++j) {
                    if (i === j) continue;
                    const d = (Y[i][0] - Y[j][0]) ** 2 + (Y[i][1] - Y[j][1]) ** 2;
                    if (d < bestD) {
                        bestD = d;
                        best = j;
                    }
                }
                if (labels[best] === labels[i]) ++agree;
            }
            // Chance is ~1/3; a collapsed or diverged embedding scores near that.
            expect(agree / X.length).toBeGreaterThan(0.9);
        });
    });
});
