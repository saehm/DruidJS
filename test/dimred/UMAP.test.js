import { describe, expect, it } from "vitest";
import { UMAP } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

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
