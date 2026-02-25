import { describe, expect, it } from "vitest";
import { TriMap } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("TriMap", () => {
    it("should complete within timeout", { timeout: 30000 }, () => {
        const data = generateTestData(20, 5);
        const trimap = new TriMap(data, { d: 2, c: 2, seed: 42 });
        const result = trimap.transform(50);

        expect(result).toHaveLength(20);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 30000 }, () => {
        const data = generateTestData(15, 4);
        const trimap = new TriMap(data, { d: 2, c: 2, seed: 42 });
        const result = trimap.transform(30);
        expectValidValues(result, "TriMap");
    });
    it("should work with static transform", () => {
        const data = generateTestData(60, 4);
        const result = TriMap.transform(data, { d: 2 });
        expect(result).toHaveLength(60);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with static generator", () => {
        const data = generateTestData(60, 4);
        const gen = TriMap.generator(data, { d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(60);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(60, 4);
        const result = await TriMap.transform_async(data, { d: 2 });
        expect(result).toHaveLength(60);
    });
});
