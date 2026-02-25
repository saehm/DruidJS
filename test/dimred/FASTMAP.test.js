import { describe, expect, it } from "vitest";
import { FASTMAP } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("FASTMAP", () => {
    it("should reduce dimensionality to d=2 by default", () => {
        const data = generateTestData(15, 4);
        const fastmap = new FASTMAP(data, { seed: 42 });
        const result = fastmap.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with generator", () => {
        const data = generateTestData(10, 5);
        const fastmap = new FASTMAP(data);
        const gen = fastmap.generator();
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(10);
    });

    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = FASTMAP.transform(data);
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = FASTMAP.generator(data);
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await FASTMAP.transform_async(data);
        expect(result).toHaveLength(10);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 4);
        const fastmap = new FASTMAP(data, { d: 2, seed: 42 });
        const result = fastmap.transform();
        expectValidValues(result, "FASTMAP");
    });
});
