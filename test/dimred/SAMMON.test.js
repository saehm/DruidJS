import { describe, expect, it } from "vitest";
import { SAMMON } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("SAMMON", () => {
    it("should complete within timeout", { timeout: 20000 }, () => {
        const data = generateTestData(15, 4);
        const sammon = new SAMMON(data, { d: 2, init_DR: "random", seed: 42 });
        const result = sammon.transform(50);

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 20000 }, () => {
        const data = generateTestData(12, 4);
        const sammon = new SAMMON(data, { d: 2, init_DR: "random", seed: 42 });
        const result = sammon.transform(30);
        expectValidValues(result, "SAMMON");
    });

    it("generator should yield intermediate results", { timeout: 20000 }, () => {
        const data = generateTestData(12, 4);
        const sammon = new SAMMON(data, { d: 2, init_DR: "random", seed: 42 });
        const gen = sammon.generator(10);

        let count = 0;
        for (const result of gen) {
            count++;
            expect(result).toHaveLength(12);
            expect(result[0]).toHaveLength(2);
        }
        expect(count).toBe(10);
    });
    it("should use default parameters if none provided", () => {
        const data = generateTestData(15, 4);
        const sammon = new SAMMON(data);
        expect(sammon.parameter("magic")).toBe(0.1);
    });

    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = SAMMON.transform(data, { d: 2 });
        expect(result).toHaveLength(10);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = SAMMON.generator(data, { d: 2 });
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await SAMMON.transform_async(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should handle PCA initialization", () => {
        const data = generateTestData(10, 5);
        const sammon = new SAMMON(data, { init_DR: "PCA", d: 2 });
        const result = sammon.transform(10);
        expect(result).toHaveLength(10);
    });

    it("should handle MDS initialization", () => {
        const data = generateTestData(10, 5);
        const sammon = new SAMMON(data, { init_DR: "MDS", d: 2 });
        const result = sammon.transform(10);
        expect(result).toHaveLength(10);
    });

    it("should throw error for invalid init_DR", () => {
        const data = generateTestData(10, 5);
        const sammon = new SAMMON(data, { init_DR: "invalid" });
        expect(() => sammon.init()).toThrow();
    });
});
