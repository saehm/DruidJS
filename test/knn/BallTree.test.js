import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { BallTree } from "../../src/knn/index.js";
import { euclidean, manhattan } from "../../src/metrics/index.js";

describe("BallTree", () => {
    it("should find nearest neighbors in 2D", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
            [10, 10],
        ];
        const tree = new BallTree(points, { metric: euclidean });

        const neighbors = tree.search([0, 0], 2);
        expect(neighbors).toHaveLength(2);
        expect(neighbors[0].index).toBe(0);
    });

    it("should work with IRIS dataset", () => {
        const points = mistle.IRIS.values;
        const tree = new BallTree(points, { metric: euclidean });

        const neighbors = tree.search(points[0], 3);
        expect(neighbors).toHaveLength(3);
        expect(neighbors[0].index).toBe(0);
    });

    it("should handle large dataset (YEAST)", () => {
        const points = mistle.YEAST.values;
        const tree = new BallTree(points, { metric: euclidean });
        const neighbors = tree.search(points[0], 10);

        expect(neighbors).toHaveLength(10);
        expect(neighbors[0].index).toBe(0);
    });

    it("should work with different metrics", () => {
        const points = mistle.WINE.values;
        const tree = new BallTree(points, { metric: manhattan });
        const neighbors = tree.search(points[0], 5);

        expect(neighbors).toHaveLength(5);
        expect(neighbors[0].distance).toBe(0);
    });

    it("should work with search_by_index", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
        ];
        const tree = new BallTree(points, { metric: euclidean });
        const neighbors = tree.search_by_index(0, 2);
        expect(neighbors).toHaveLength(2);
        expect(neighbors[0].index).toBe(0);
    });
});
