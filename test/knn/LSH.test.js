import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { LSH } from "../../src/knn/index.js";
import { euclidean } from "../../src/metrics/index.js";

describe("LSH", () => {
    it("should find nearest neighbors in 2D", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
            [10, 10],
        ];
        const lsh = new LSH(points, { metric: euclidean });

        const neighbors = lsh.search([0, 0], 2);
        expect(neighbors).toHaveLength(2);
        expect(neighbors[0].index).toBe(0);
    });

    it("should work with IRIS dataset", () => {
        const points = mistle.IRIS.values;
        const lsh = new LSH(points, { metric: euclidean });

        const neighbors = lsh.search(points[0], 3);
        expect(neighbors).toHaveLength(3);
        expect(neighbors[0].index).toBe(0);
    });

    it("should handle large dataset (YEAST)", () => {
        const points = mistle.YEAST.values;
        const lsh = new LSH(points, { metric: euclidean });
        const neighbors = lsh.search(points[0], 10);

        expect(neighbors).toHaveLength(10);
    });

    it("should work with search_by_index", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
        ];
        const lsh = new LSH(points, { metric: euclidean });
        const neighbors = lsh.search_by_index(0, 2);
        expect(neighbors).toHaveLength(2);
        expect(neighbors[0].index).toBe(0);
    });

    it("should fallback to linear search when candidates are insufficient", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
        ];
        // Use 1 hash table to increase chance of insufficient candidates for large k
        const lsh = new LSH(points, { metric: euclidean, numHashTables: 1 });
        const neighbors = lsh.search([0, 0], 5);
        expect(neighbors.length).toBeLessThanOrEqual(3);
    });
});
