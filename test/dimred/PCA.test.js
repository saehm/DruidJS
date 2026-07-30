import { describe, expect, it } from "vitest";
import { PCA } from "../../src/dimred/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { generateTestData } from "../utils/data-generators.js";

describe("PCA", () => {
    it("should reduce dimensionality to d=2 by default", () => {
        const data = generateTestData();
        const pca = new PCA(data, { seed: 42 });
        const result = pca.transform();

        expect(result).toHaveLength(10);
        expect(result[0]).toHaveLength(2);
    });

    it("should reduce dimensionality to specified d", () => {
        const data = generateTestData();
        const pca = new PCA(data, { d: 3, seed: 42 });
        const result = pca.transform();

        expect(result).toHaveLength(10);
        expect(result[0]).toHaveLength(3);
    });

    it("should work with Matrix input", () => {
        const data = Matrix.from(generateTestData());
        const pca = new PCA(data, { d: 2, seed: 42 });
        const result = pca.transform();

        expect(result.shape[0]).toBe(10);
        expect(result.shape[1]).toBe(2);
    });

    it("should return principal components", () => {
        const data = generateTestData();
        const pca = new PCA(data, { d: 2, seed: 42 });
        const V = pca.principal_components();

        expect(V.shape[0]).toBe(4); // original dimensions
        expect(V.shape[1]).toBe(2); // target dimensions
    });

    it("should produce consistent results with same seed", () => {
        const data = generateTestData();
        const pca1 = new PCA(data, { d: 2, seed: 42 });
        const pca2 = new PCA(data, { d: 2, seed: 42 });

        const result1 = pca1.transform();
        const result2 = pca2.transform();

        expect(result1).toEqual(result2);
    });
    it("should work with static transform", () => {
        const data = generateTestData();
        const result = PCA.transform(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData();
        const gen = PCA.generator(data, { d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData();
        const result = await PCA.transform_async(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should cache principal components", () => {
        const data = generateTestData();
        const pca = new PCA(data, { d: 2 });
        const V1 = pca.principal_components();
        const V2 = pca.principal_components();
        expect(V1).toBe(V2);
    });

    it("should keep the cache across a transform", () => {
        const data = generateTestData();
        const pca = new PCA(data, { d: 2, seed: 42 });
        const V1 = pca.principal_components();
        pca.transform();
        expect(pca.principal_components()).toBe(V1);
    });

    it("should drop the cached components when d changes", () => {
        const data = generateTestData();
        const pca = new PCA(data, { d: 3, seed: 42 });
        expect(pca.principal_components().shape[1]).toBe(3);
        expect(/** @type {number[][]} */ (pca.transform())[0]).toHaveLength(3);

        // Without invalidation the components computed for d=3 are returned again and the
        // projection silently keeps the old dimensionality.
        pca.parameter("d", 2);
        expect(pca.principal_components().shape[1]).toBe(2);
        expect(/** @type {number[][]} */ (pca.transform())[0]).toHaveLength(2);
    });

    it("should recompute the projection when d changes", () => {
        const data = generateTestData();
        const changed = new PCA(data, { d: 3, seed: 42 });
        changed.transform();
        changed.parameter("d", 2);
        const recomputed = /** @type {number[][]} */ (changed.transform());

        const direct = /** @type {number[][]} */ (new PCA(data, { d: 2, seed: 42 }).transform());

        // Compared up to sign, and only to ~1e-5: an eigenvector's orientation is arbitrary, and
        // the two instances reach `simultaneous_poweriteration` with their randomizer at different
        // states — so the axes can come back mirrored, and the iteration stops at its convergence
        // tolerance rather than at the exact eigenvector. The coordinates are what must agree.
        expect(recomputed).toHaveLength(direct.length);
        for (let i = 0; i < direct.length; ++i) {
            for (let j = 0; j < 2; ++j) {
                expect(Math.abs(recomputed[i][j])).toBeCloseTo(Math.abs(direct[i][j]), 4);
            }
        }
    });

    it("should work with static principal_components", () => {
        const data = generateTestData();
        const V = PCA.principal_components(data, { d: 2 });
        expect(V.shape[0]).toBe(4);
        expect(V.shape[1]).toBe(2);
    });
});
