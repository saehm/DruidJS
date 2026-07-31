import { describe, expect, it } from "vitest";
import { LLE } from "../../src/dimred/index.js";
import { HNSW, NaiveKNN, spatial_tree } from "../../src/knn/index.js";
import { euclidean } from "../../src/metrics/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("LLE", () => {
    it("should reduce dimensionality to specified d", () => {
        const data = generateTestData(15, 4);
        const lle = new LLE(data, { neighbors: 4, d: 2, seed: 42 });
        const result = lle.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 4);
        const lle = new LLE(data, { neighbors: 4, d: 2, seed: 42 });
        const result = lle.transform();
        expectValidValues(result, "LLE");
    });
    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = LLE.transform(data, { neighbors: 4, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 4);
        const gen = LLE.generator(data, { neighbors: 2, d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 4);
        const result = await LLE.transform_async(data, { neighbors: 2, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should use default neighbors", () => {
        const data = generateTestData(20, 5);
        const lle = new LLE(data);
        expect(lle.parameter("neighbors")).toBe(2);
    });

    it("should accept an exact knn index and match the default", () => {
        const data = generateTestData(20, 4);
        const expected = new LLE(data, { neighbors: 4, d: 2, seed: 42 }).transform();

        for (const knn of [
            spatial_tree(data, { metric: euclidean, seed: 42 }),
            new NaiveKNN(data, { metric: euclidean, seed: 42 }),
        ]) {
            const actual = new LLE(data, { neighbors: 4, d: 2, seed: 42, knn }).transform();
            expect(actual).toEqual(expected);
        }
    });

    it("should accept an approximate knn index", () => {
        const data = generateTestData(20, 4);
        const knn = new HNSW(data, { metric: euclidean, seed: 42, ef: 100 });
        const result = new LLE(data, { neighbors: 4, d: 2, seed: 42, knn }).transform();
        expect(result).toHaveLength(20);
        expectValidValues(result, "LLE");
    });

    it("should apply regularization when neighbors > dimensions", () => {
        const data = generateTestData(10, 2); // 10 points, 2 dims
        const lle = new LLE(data, { neighbors: 5, d: 1 });
        const result = lle.transform();
        expect(result).toHaveLength(10);
    });
});
