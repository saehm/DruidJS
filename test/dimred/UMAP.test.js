import { describe, expect, it } from "vitest";
import { UMAP } from "../../src/dimred/index.js";
import { spatial_tree } from "../../src/knn/index.js";
import { euclidean } from "../../src/metrics/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("UMAP._smooth_knn_dist", () => {
    // Values cross-checked against umap-learn 0.5.12's smooth_knn_dist on the same input.
    const data = generateTestData(60, 3);
    const k = 10;
    const tree = spatial_tree(data, { metric: euclidean, seed: 1212 });
    const umap = new UMAP(data, { n_neighbors: k, local_connectivity: 1, d: 2, seed: 1212 });
    // @ts-ignore - private
    const { sigmas, rhos } = umap._smooth_knn_dist(tree, k);

    // The row UMAP works on: the point itself at distance 0, then its k-1 nearest neighbors.
    const rowOf = (i) => [0, ...tree.search_by_index(i, k - 1).map((n) => n.distance)];

    it("should set rho to the distance to the nearest neighbor", () => {
        for (let i = 0; i < data.length; ++i) {
            const nearest = rowOf(i).filter((d) => d > 0);
            expect(rhos[i]).toBeCloseTo(Math.min(...nearest), 10);
        }
    });

    it("should converge sigma rather than clamping to the MIN_K_DIST_SCALE floor", () => {
        // The binary search used to never converge, leaving every sigma pinned at 1e-3 * mean_d.
        for (let i = 0; i < data.length; ++i) {
            const row = rowOf(i);
            const floor = 1e-3 * (row.reduce((a, b) => a + b, 0) / row.length);
            expect(sigmas[i]).toBeGreaterThan(floor * 10);
            expect(Number.isFinite(sigmas[i])).toBe(true);
        }
    });

    it("should solve sum_j exp(-(d_j - rho)/sigma) = log2(k) over the neighbors", () => {
        const target = Math.log2(k);
        for (let i = 0; i < data.length; ++i) {
            const row = rowOf(i);
            let psum = 0;
            // Starts at 1: the self entry is not part of the fuzzy set.
            for (let j = 1; j < row.length; ++j) {
                const d = row[j] - rhos[i];
                psum += d > 0 ? Math.exp(-(d / sigmas[i])) : 1;
            }
            expect(psum).toBeCloseTo(target, 4);
        }
    });
});

describe("UMAP", () => {
    it("should reduce dimensionality to d=2 by default", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const umap = new UMAP(data, { n_neighbors: 4, seed: 42 });
        const result = umap.transform(50);

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should complete within timeout", { timeout: 30000 }, () => {
        const data = generateTestData(20, 5);
        const umap = new UMAP(data, { n_neighbors: 5, d: 2, seed: 42 });
        const result = umap.transform(50);

        expect(result).toHaveLength(20);
        expect(result[0]).toHaveLength(2);
    });

    it("generator should yield intermediate results", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const umap = new UMAP(data, { n_neighbors: 5, d: 2, seed: 42 });
        const gen = umap.generator(10);

        let count = 0;
        for (const result of gen) {
            count++;
            expect(result).toHaveLength(15);
            expect(result[0]).toHaveLength(2);
        }
        expect(count).toBe(10);
    });

    it("should produce valid values", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const umap = new UMAP(data, { n_neighbors: 4, d: 2, seed: 42 });
        const result = umap.transform(50);
        expectValidValues(result, "UMAP");
    });
    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = UMAP.transform(data, { n_neighbors: 4 });
        expect(result).toHaveLength(10);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = UMAP.generator(data, { n_neighbors: 4 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await UMAP.transform_async(data, { n_neighbors: 4 });
        expect(result).toHaveLength(10);
    });
});
