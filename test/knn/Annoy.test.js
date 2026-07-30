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
        // search_by_index excludes the queried element itself.
        expect(neighbors.map((n) => n.index)).toEqual([1, 2]);
    });

    it("should provide search_index alias", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const annoy = new Annoy(points, { seed: 42 });
        expect(annoy.search_index(0, 2)).toEqual(annoy.search_by_index(0, 2));
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

    describe("built empty and filled with add()", () => {
        const points = mistle.IRIS.values;

        it("should report indices into the added elements, not shifted by the dummy", () => {
            const annoy = new Annoy([], { metric: euclidean, seed: 42 });
            annoy.add(points);

            // The empty constructor needs an element to get past KNN's non-empty check. If that
            // placeholder survives into `_elements` it sits at index 0 and shifts everything.
            const neighbors = annoy.search(points[7], 5);
            expect(neighbors[0].index).toBe(7);
            expect(neighbors[0].distance).toBe(0);
            for (const { index, element } of neighbors) {
                expect(Array.from(element)).toEqual(Array.from(points[index]));
            }
        });

        it("should match an index constructed with the same points directly", () => {
            const direct = new Annoy(points, { metric: euclidean, seed: 42 });
            const incremental = new Annoy([], { metric: euclidean, seed: 42 });
            incremental.add(points);

            expect(incremental.search(points[3], 5).map((n) => n.index)).toEqual(
                direct.search(points[3], 5).map((n) => n.index),
            );
        });

        it("should not admit a phantom point at the origin", () => {
            const far = [
                [100, 100],
                [101, 101],
                [102, 102],
            ];
            const annoy = new Annoy([], { metric: euclidean, seed: 42 });
            annoy.add(far);

            // A dummy `new Float64Array([0])` left in place would be the nearest thing to [0, 0].
            const neighbors = annoy.search([0, 0], 3);
            expect(neighbors).toHaveLength(3);
            for (const { element } of neighbors) {
                expect(Array.from(element)).not.toEqual([0]);
            }
        });

        it("should keep indices contiguous across several add() calls", () => {
            const annoy = new Annoy([], { metric: euclidean, seed: 42 });
            annoy.add(points.slice(0, 40));
            annoy.add(points.slice(40));

            const neighbors = annoy.search(points[95], 3);
            expect(neighbors[0].index).toBe(95);
            expect(neighbors[0].distance).toBe(0);
        });

        it("should search an index that was never filled", () => {
            const annoy = new Annoy([], { metric: euclidean, seed: 42 });
            expect(annoy.search([0, 0], 3)).toEqual([]);
        });
    });
});
