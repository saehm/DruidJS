import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { Annoy } from "../../src/knn/index.js";
import { euclidean, manhattan } from "../../src/metrics/index.js";

describe("Annoy", () => {
    it("should find nearest neighbors in 2D", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
            [10, 10],
        ];
        const annoy = new Annoy(points, { metric: euclidean });

        const neighbors = annoy.search([0, 0], 2);
        expect(neighbors).toHaveLength(2);
        expect(neighbors[0].index).toBe(0);
    });

    it("should work with IRIS dataset", () => {
        const points = mistle.IRIS.values;
        const annoy = new Annoy(points, { metric: euclidean });

        const neighbors = annoy.search(points[0], 3);
        expect(neighbors).toHaveLength(3);
        expect(neighbors[0].index).toBe(0);
    });

    it("should handle large dataset (YEAST)", () => {
        const points = mistle.YEAST.values;
        const annoy = new Annoy(points, { metric: euclidean });
        const neighbors = annoy.search(points[0], 10);

        expect(neighbors).toHaveLength(10);
        expect(neighbors[0].index).toBe(0);
    });

    it("should work with different metrics", () => {
        const points = mistle.WINE.values;
        const annoy = new Annoy(points, { metric: manhattan });
        const neighbors = annoy.search(points[0], 5);

        expect(neighbors).toHaveLength(5);
        expect(neighbors[0].distance).toBe(0);
    });

    it("should work with search_by_index", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
        ];
        const annoy = new Annoy(points, { metric: euclidean });
        const neighbors = annoy.search_by_index(0, 2);
        expect(neighbors).toHaveLength(2);
        expect(neighbors[0].index).toBe(0);
    });

    it("should provide search_index alias", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const annoy = new Annoy(points, { seed: 42 });
        const neighbors = annoy.search_index(0, 2);
        expect(neighbors).toHaveLength(2);
    });

    it("should have correct num_trees and num_nodes", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const annoy = new Annoy(points, { n_trees: 5, seed: 42 });
        expect(annoy.num_trees).toBeGreaterThanOrEqual(5);
        expect(annoy.num_nodes).toBeGreaterThan(0);
    });

    it("should fall back to linear search if k is large", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
            [3, 3],
            [4, 4],
        ];
        const annoy = new Annoy(points, { n_trees: 1, maxPointsPerLeaf: 2, seed: 42 });
        // Requesting more neighbors than exist total should definitely trigger fallback
        const results = annoy.search([0.5, 0.5], 10);
        expect(results.length).toBeLessThanOrEqual(5);
    });
});
