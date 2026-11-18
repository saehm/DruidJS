import { describe, expect, it } from "vitest";
import { FASTMAP } from "../../src/dimred/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { cosine, manhattan } from "../../src/metrics/index.js";
import { Randomizer } from "../../src/util/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("FASTMAP", () => {
    it("should reduce dimensionality to d=2 by default", () => {
        const data = generateTestData(15, 4);
        const fastmap = new FASTMAP(data, { seed: 42 });
        const result = fastmap.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with generator", () => {
        const data = generateTestData(10, 5);
        const fastmap = new FASTMAP(data);
        const gen = fastmap.generator();
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(10);
    });

    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = FASTMAP.transform(data);
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = FASTMAP.generator(data);
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await FASTMAP.transform_async(data);
        expect(result).toHaveLength(10);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 4);
        const fastmap = new FASTMAP(data, { d: 2, seed: 42 });
        const result = fastmap.transform();
        expectValidValues(result, "FASTMAP");
    });

    describe("distances", () => {
        const r = new Randomizer(2);
        const X = new Matrix(200, 6, () => r.random);

        it("is a linear projection for Euclidean distances", () => {
            // Y = [X, 1] B exactly: every axis is the projection onto a line through two pivots
            const Y = /** @type {Matrix} */ (new FASTMAP(X, { d: 5 }).transform());
            const A = new Matrix(200, 7, (i, j) => (j < 6 ? X.entry(i, j) : 1));
            const B = Matrix.solve(A.transDot(A), A.transDot(Y));
            const R = A.dot(B).sub(Y);
            expect(Math.max(...R.values.map(Math.abs))).toBeLessThan(1e-8);
        });

        // The leftover squared distance d² − (y_a − y_b)² can go negative for non-Euclidean
        // distances; its square root used to be NaN, and the NaNs then broke the pivot search.
        it.each([
            ["manhattan", manhattan],
            ["cosine", cosine],
        ])("stays finite with %s distances on later axes", (_, metric) => {
            const Y = /** @type {Matrix} */ (new FASTMAP(X, { metric, d: 4 }).transform());
            expect(Y.values.every(Number.isFinite)).toBe(true);
        });
    });
});
