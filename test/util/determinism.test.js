import { describe, expect, it } from "vitest";
import { ISOMAP } from "../../src/dimred/ISOMAP.js";
import { TriMap } from "../../src/dimred/TriMap.js";
import { UMAP } from "../../src/dimred/UMAP.js";
import { BallTree } from "../../src/knn/BallTree.js";
import { KDTree } from "../../src/knn/KDTree.js";
import { spatial_tree } from "../../src/knn/spatial_tree.js";
import { cosine, euclidean, manhattan } from "../../src/metrics/index.js";
import { quickselect, quickselectByAxis } from "../../src/util/quickselect.js";
import { Randomizer } from "../../src/util/randomizer.js";

/**
 * Everything in the library derives its randomness from a seeded `Randomizer`, so the same seed and
 * the same input must always give the same output. `quickselect` picking pivots with `Math.random`
 * silently broke that: the pivot only changes the outcome when keys tie on the split axis, which is
 * exactly what happens with integer, one-hot, or quantized features — and never with the continuous
 * random data the rest of the suite uses.
 */

/** Coordinates with heavy ties on every axis. */
function tied_points(n = 200) {
    return Array.from({ length: n }, (_, i) => [i % 5, i % 3, i % 7]);
}

/** A fresh randomizer per call, so each assertion is independent of the ones before it. */
const seeded = () => new Randomizer(1212);

describe("quickselect is deterministic", () => {
    it("returns the same partition for the same input every time", () => {
        const input = [9, 3, 7, 1, 5, 2, 8, 4, 6, 3, 7, 1];
        const results = new Set();
        for (let run = 0; run < 16; ++run) {
            const copy = [...input];
            quickselect(copy, seeded(), 5);
            results.add(JSON.stringify(copy));
        }
        expect(results.size).toBe(1);
    });

    it("draws pivots from the randomizer, never from Math.random", () => {
        // Patching Math.random would change the outcome if quickselect still reached for it.
        const original = Math.random;
        Math.random = () => {
            throw new Error("quickselect must not use Math.random");
        };
        try {
            const arr = Array.from({ length: 300 }, (_, i) => (i * 37) % 91);
            expect(() => quickselect(arr, seeded(), 150)).not.toThrow();
        } finally {
            Math.random = original;
        }
    });

    it("selects the correct rank even when every key is equal", () => {
        const arr = new Array(500).fill(7);
        expect(quickselect(arr, seeded(), 250)).toBe(7);
    });

    it("selects the correct rank with many duplicates", () => {
        const arr = Array.from({ length: 300 }, (_, i) => i % 10);
        const sorted = [...arr].sort((a, b) => a - b);
        for (const k of [0, 42, 150, 299]) {
            expect(quickselect([...arr], seeded(), k)).toBe(sorted[k]);
        }
    });

    it("partitions by axis around the true k-th element", () => {
        const points = tied_points(101);
        const arr = points.map((element, index) => ({ index, element }));
        const k = 50;
        quickselectByAxis(arr, seeded(), k, 0);
        const pivot = arr[k].element[0];
        for (let i = 0; i < k; ++i) expect(arr[i].element[0]).toBeLessThanOrEqual(pivot);
        for (let i = k + 1; i < arr.length; ++i) expect(arr[i].element[0]).toBeGreaterThanOrEqual(pivot);
    });
});

describe("spatial trees are reproducible on tied data", () => {
    for (const [name, Tree] of [
        ["KDTree", KDTree],
        ["BallTree", BallTree],
    ]) {
        it(`${name} returns identical neighbors across rebuilds`, () => {
            const results = new Set();
            for (let run = 0; run < 8; ++run) {
                const tree = new Tree(tied_points(), { metric: euclidean, seed: 1212 });
                results.add(JSON.stringify(tree.search([2, 1, 3], 5).map((d) => d.index)));
            }
            expect(results.size).toBe(1);
        });
    }

    it("defaults to the euclidean metric when only a seed is given", () => {
        expect(() => new KDTree(tied_points(), { seed: 1212 }).search([2, 1, 3], 3)).not.toThrow();
        expect(() => new BallTree(tied_points(), { seed: 1212 }).search([2, 1, 3], 3)).not.toThrow();
    });
});

describe("spatial_tree picks an index that is correct for the metric", () => {
    it("uses a KDTree for L_p metrics", () => {
        expect(spatial_tree(tied_points(), { metric: euclidean, seed: 1212 })).toBeInstanceOf(KDTree);
        expect(spatial_tree(tied_points(), { metric: manhattan, seed: 1212 })).toBeInstanceOf(KDTree);
    });

    it("falls back to a BallTree for metrics a KDTree cannot prune soundly", () => {
        expect(spatial_tree(tied_points(), { metric: cosine, seed: 1212 })).toBeInstanceOf(BallTree);
    });

    it("agrees with brute force for a non-L_p metric", () => {
        const points = Array.from({ length: 120 }, (_, i) => [((i * 7) % 13) + 1, ((i * 5) % 11) + 1, (i % 9) + 1]);
        const query = [4, 6, 3];
        const tree = spatial_tree(points, { metric: cosine, seed: 1212 });
        const found = tree
            .search(query, 5)
            .map((d) => d.distance)
            .sort((a, b) => a - b);
        const expected = points
            .map((p) => cosine(p, query))
            .sort((a, b) => a - b)
            .slice(0, 5);
        for (let i = 0; i < 5; ++i) expect(found[i]).toBeCloseTo(expected[i], 12);
    });
});

describe("dimensionality reduction is reproducible on tied data", () => {
    const X = Array.from({ length: 120 }, (_, i) => [i % 6, i % 4, i % 5, i % 3]);

    for (const [name, run] of [
        ["UMAP", () => new UMAP(X, { seed: 42, d: 2 }).transform(15)],
        ["ISOMAP", () => new ISOMAP(X, { seed: 42, d: 2, neighbors: 8 }).transform()],
        ["TriMap", () => new TriMap(X, { seed: 42, d: 2 }).transform(15)],
    ]) {
        it(`${name} gives the same embedding for the same seed`, () => {
            const results = new Set();
            for (let run_index = 0; run_index < 3; ++run_index) {
                results.add(JSON.stringify(run()));
            }
            expect(results.size).toBe(1);
        });
    }
});
