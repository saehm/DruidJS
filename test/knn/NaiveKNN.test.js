import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { NaiveKNN } from "../../src/knn/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { euclidean } from "../../src/metrics/index.js";
import { generateClusteredData, generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("NaiveKNN", () => {
    it("should find nearest neighbors in 2D", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
            [10, 10],
        ];
        const naive = new NaiveKNN(points, { metric: euclidean });

        const neighbors = naive.search([0, 0], 2);
        expect(neighbors).toHaveLength(2);
        expect(neighbors[0].index).toBe(0);
    });

    it("should work with IRIS dataset", () => {
        const points = mistle.IRIS.values;
        const naive = new NaiveKNN(points, { metric: euclidean });

        const neighbors = naive.search(points[0], 3);
        expect(neighbors).toHaveLength(3);
        expect(neighbors[0].index).toBe(0);
    });

    it("should handle large dataset (YEAST)", () => {
        const points = mistle.YEAST.values;
        const naive = new NaiveKNN(points, { metric: euclidean });
        const neighbors = naive.search(points[0], 10);

        expect(neighbors).toHaveLength(10);
        expect(neighbors[0].index).toBe(0);
    });

    it("should produce valid values", () => {
        const data = generateClusteredData();
        const knn = new NaiveKNN(data, { neighbors: 5, seed: 42 });
        const result = knn.search(data[0], 5);
        expectValidValues(result, "NaiveKNN");
    });

    it("should work with precomputed metric", () => {
        const data = generateTestData(10, 4);
        const dists = Matrix.from(generateTestData(10, 10)); // dummy dist matrix
        const knn = new NaiveKNN(dists, { metric: "precomputed" });
        const result = knn.search_by_index(0, 5);
        expect(result).toHaveLength(5);
    });

    it("should work with search_by_index without precomputed metric", () => {
        const data = generateTestData(10, 4);
        const knn = new NaiveKNN(data, { metric: euclidean });
        const result = knn.search_by_index(0, 5);
        expect(result).toHaveLength(5);
    });

    it("should throw error if search is called with precomputed metric", () => {
        const data = generateTestData(10, 4);
        const dists = Matrix.from(generateTestData(10, 10)); // dummy dist matrix
        const knn = new NaiveKNN(dists, { metric: "precomputed" });
        expect(() => knn.search(data[0], 5)).toThrow(
            "Search by query element is only possible when not using a precomputed distance matrix!",
        );
    });

    it("should detect typed arrays", () => {
        const data = [new Float64Array([1, 2]), new Float64Array([3, 4])];
        const knn = new NaiveKNN(data);
        // @ts-ignore
        expect(knn._type).toBe("typed");
    });

    it("should throw error if elements are empty", () => {
        expect(() => new NaiveKNN([])).toThrow("Elements needs to contain at least one element!");
    });

    it("should return the k nearest in ascending order, matching a brute force scan", () => {
        const data = generateTestData(60, 3);
        const knn = new NaiveKNN(data, { metric: euclidean, seed: 42 });
        const ranked = data
            .map((p, index) => ({ index, distance: euclidean(data[7], p) }))
            .sort((a, b) => a.distance - b.distance || a.index - b.index);

        // `search` has no notion of self: an exact match is a legitimate result.
        const bySearch = knn.search(data[7], 9);
        expect(bySearch.map((r) => r.index)).toEqual(ranked.slice(0, 9).map((r) => r.index));
        bySearch.forEach((r, i) => {
            expect(r.distance).toBeCloseTo(ranked[i].distance, 12);
        });

        // `search_by_index` drops element 7 itself and fills the gap with the next nearest.
        const expectedByIndex = ranked.filter((r) => r.index !== 7).slice(0, 9);
        const byIndex = knn.search_by_index(7, 9);
        expect(byIndex.map((r) => r.index)).toEqual(expectedByIndex.map((r) => r.index));
        byIndex.forEach((r, i) => {
            expect(r.distance).toBeCloseTo(expectedByIndex[i].distance, 12);
        });
    });

    it("should clamp k to the number of elements", () => {
        const data = generateTestData(6, 2);
        const knn = new NaiveKNN(data, { metric: euclidean });
        expect(knn.search(data[0], 100)).toHaveLength(6);
        expect(knn.search(data[0], 0)).toHaveLength(0);
    });

    it("should break ties deterministically by index", () => {
        const identical = Array.from({ length: 20 }, () => [1, 1]);
        const knn = new NaiveKNN(identical, { metric: euclidean });
        expect(knn.search([0, 0], 4).map((r) => r.index)).toEqual([0, 1, 2, 3]);
        expect(knn.search([0, 0], 4).map((r) => r.index)).toEqual([0, 1, 2, 3]);
    });

    it("should not build the distance matrix unless search_by_index is used", () => {
        const data = generateTestData(20, 3);
        const knn = new NaiveKNN(data, { metric: euclidean });
        // @ts-ignore - private
        expect(knn._D).toBeNull();
        knn.search(data[0], 5);
        // @ts-ignore - private
        expect(knn._D).toBeNull();
        knn.search_by_index(0, 5);
        // @ts-ignore - private
        expect(knn._D).toBeInstanceOf(Matrix);
    });
});
