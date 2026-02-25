import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { NNDescent } from "../../src/knn/index.js";
import { euclidean } from "../../src/metrics/index.js";

describe("NNDescent", () => {
    it("should find nearest neighbors in 2D", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
            [10, 10],
        ];
        const nnd = new NNDescent(points, {
            metric: euclidean,
            samples: 2,
            seed: 42,
        });

        const neighbors = nnd.search([0, 0], 2);
        expect(neighbors).toHaveLength(2);
        expect(neighbors[0].index).toBe(0);
    });

    it("should work with IRIS dataset", () => {
        const points = mistle.IRIS.values;
        const nnd = new NNDescent(points, {
            metric: euclidean,
            samples: 5,
            seed: 42,
        });

        const neighbors = nnd.search(points[0], 3);
        expect(neighbors).toHaveLength(3);
        expect(neighbors[0].index).toBe(0);
    });

    it("should handle large dataset (YEAST)", () => {
        const points = mistle.YEAST.values;
        const nnd = new NNDescent(points, {
            metric: euclidean,
            samples: 10,
            seed: 42,
        });

        const neighbors = nnd.search(points[0], 10);
        expect(neighbors).toHaveLength(10);
        expect(neighbors[0].index).toBe(0);
    });

    it("should work with search_by_index and search_index", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
        ];
        const nnd = new NNDescent(points, { samples: 2, seed: 42 });
        const neighbors1 = nnd.search_by_index(0, 2);
        expect(neighbors1).toHaveLength(2);
        expect(neighbors1[0].index).toBe(0);
    });

    it("should handle invalid search_by_index indices", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const nnd = new NNDescent(points, { samples: 2 });
        expect(nnd.search_by_index(-1)).toHaveLength(0);
        expect(nnd.search_by_index(10)).toHaveLength(0);
    });
});
