import { describe, expect, it } from "vitest";
import { PaCMAP } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("PaCMAP", () => {
    it("should reduce dimensionality to d=2 by default", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const pacmap = new PaCMAP(data, { n_neighbors: 4, seed: 42 });
        const result = pacmap.transform(50);

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should complete within timeout", { timeout: 30000 }, () => {
        const data = generateTestData(20, 5);
        const pacmap = new PaCMAP(data, { n_neighbors: 5, d: 2, seed: 42 });
        const result = pacmap.transform(50);

        expect(result).toHaveLength(20);
        expect(result[0]).toHaveLength(2);
    });

    it("generator should yield intermediate results", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const pacmap = new PaCMAP(data, { n_neighbors: 5, d: 2, seed: 42 });
        const gen = pacmap.generator(10);

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
        const pacmap = new PaCMAP(data, { n_neighbors: 4, d: 2, seed: 42 });
        const result = pacmap.transform(50);
        expectValidValues(result, "PaCMAP");
    });

    it("should work with static transform", { timeout: 30000 }, () => {
        const data = generateTestData(12, 5);
        const result = PaCMAP.transform(data, { n_neighbors: 4 });
        expect(result).toHaveLength(12);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with static generator", { timeout: 30000 }, () => {
        const data = generateTestData(12, 5);
        const gen = PaCMAP.generator(data, { n_neighbors: 4 });
        const result = gen.next().value;
        expect(result).toHaveLength(12);
    });

    it("should work with static transform_async", { timeout: 30000 }, async () => {
        const data = generateTestData(12, 5);
        const result = await PaCMAP.transform_async(data, { n_neighbors: 4 });
        expect(result).toHaveLength(12);
    });

    it("should respect d=3 parameter", { timeout: 30000 }, () => {
        const data = generateTestData(15, 6);
        const pacmap = new PaCMAP(data, { n_neighbors: 4, d: 3, seed: 42 });
        const result = pacmap.transform(30);

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(3);
    });

    it("should throw if n_neighbors >= N", () => {
        const data = generateTestData(10, 4);
        expect(() => new PaCMAP(data, { n_neighbors: 10 })).toThrow();
    });
});
