import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { LTSA } from "../../src/dimred/index.js";
import { HNSW, NaiveKNN, spatial_tree } from "../../src/knn/index.js";
import { euclidean } from "../../src/metrics/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { bestRankCorrelation, expectValidValues } from "../utils/helpers.js";

describe("LTSA", () => {
    it("should unroll the Swiss roll like scikit-learn", () => {
        // scikit-learn 1.7.2 scores 1.00 / 1.00 on roll / height here. Before the fix (wrong
        // centering, wrong constant term, largest instead of smallest eigenvectors) DruidJS scored
        // 0.72 / 0.37.
        const X = mistle.swissroll({ N: 800, noise: 0.05, seed: 1212 }).values;
        const Y = new LTSA(X, { d: 2, neighbors: 10 }).transform();
        expectValidValues(Y, "LTSA");
        expect(bestRankCorrelation(Y, X.map((_, i) => i))).toBeGreaterThan(0.99);
        expect(bestRankCorrelation(Y, X.map((row) => row[2]))).toBeGreaterThan(0.99);
    });
    it("should reduce dimensionality to specified d", () => {
        const data = generateTestData(15, 5);
        const ltsa = new LTSA(data, { neighbors: 4, d: 2, seed: 42 });
        const result = ltsa.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 5);
        const ltsa = new LTSA(data, { neighbors: 4, d: 2, seed: 42 });
        const result = ltsa.transform();
        expectValidValues(result, "LTSA");
    });
    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = LTSA.transform(data, { neighbors: 4, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = LTSA.generator(data, { neighbors: 4, d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await LTSA.transform_async(data, { neighbors: 4, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should use default neighbors", () => {
        const data = generateTestData(20, 5);
        const ltsa = new LTSA(data);
        expect(ltsa.parameter("neighbors")).toBe(2);
    });

    it("should accept an exact knn index and match the default", () => {
        const data = generateTestData(20, 5);
        const expected = new LTSA(data, { neighbors: 4, d: 2, seed: 42 }).transform();

        for (const knn of [
            spatial_tree(data, { metric: euclidean, seed: 42 }),
            new NaiveKNN(data, { metric: euclidean, seed: 42 }),
        ]) {
            const actual = new LTSA(data, { neighbors: 4, d: 2, seed: 42, knn }).transform();
            expect(actual).toEqual(expected);
        }
    });

    it("should accept an approximate knn index", () => {
        const data = generateTestData(20, 5);
        const knn = new HNSW(data, { metric: euclidean, seed: 42, ef: 100 });
        const result = new LTSA(data, { neighbors: 4, d: 2, seed: 42, knn }).transform();
        expect(result).toHaveLength(20);
        expectValidValues(result, "LTSA");
    });

    it("should throw error if d >= D", () => {
        const data = generateTestData(10, 2);
        expect(() => new LTSA(data, { d: 2 })).toThrow("Dimensionality of X (D = 2) must be greater than");
    });
});
