import { describe, expect, it } from "vitest";
import { TopoMap } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("TopoMap", () => {
    it("should complete within timeout", { timeout: 10000 }, () => {
        const data = generateTestData(15, 4);
        const topomap = new TopoMap(data, { seed: 42 });
        const result = topomap.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with generator", () => {
        const data = generateTestData(10, 5);
        const topomap = new TopoMap(data);
        const gen = topomap.generator();
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(10);
    });

    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = TopoMap.transform(data);
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = TopoMap.generator(data);
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await TopoMap.transform_async(data);
        expect(result).toHaveLength(10);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(12, 4);
        const topomap = new TopoMap(data, { seed: 42 });
        const result = topomap.transform();
        expectValidValues(result, "TopoMap");
    });
});
