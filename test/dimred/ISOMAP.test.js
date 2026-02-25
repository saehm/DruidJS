import { describe, expect, it } from "vitest";
import { ISOMAP } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("ISOMAP", () => {
    it("should reduce dimensionality to d=2 by default", () => {
        const data = generateTestData(15, 4);
        const isomap = new ISOMAP(data, { neighbors: 4, seed: 42 });
        const result = isomap.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 4);
        const isomap = new ISOMAP(data, { neighbors: 4, d: 2, seed: 42 });
        const result = isomap.transform();
        expectValidValues(result, "ISOMAP");
    });
    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = ISOMAP.transform(data, { neighbors: 4, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 4);
        const gen = ISOMAP.generator(data, { neighbors: 2, d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 4);
        const result = await ISOMAP.transform_async(data, { neighbors: 2, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should use default neighbors if not provided", () => {
        const data = generateTestData(20, 5);
        const isomap = new ISOMAP(data);
        expect(isomap.parameter("neighbors")).toBe(2);
    });

    it("should work with SMACOF projection", () => {
        const data = generateTestData(15, 4);
        const isomap = new ISOMAP(data, { neighbors: 4, seed: 42, project: "SMACOF" });
        const result = isomap.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });
});
