import { describe, expect, it } from "vitest";
import { LTSA } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("LTSA", () => {
    it("should reduce dimensionality to specified d", () => {
        const data = generateTestData(15, 5);
        const ltsa = new LTSA(data, { neighbors: 4, d: 2, seed: 42 });
        const result = ltsa.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 5);
        const ltsa = new LTSA(data, { neighbors: 4, d: 2, seed: 42 });
        const result = ltsa.transform();
        expectValidValues(result, "LTSA");
    });
    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = LTSA.transform(data, { neighbors: 4, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = LTSA.generator(data, { neighbors: 4, d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await LTSA.transform_async(data, { neighbors: 4, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should use default neighbors", () => {
        const data = generateTestData(20, 5);
        const ltsa = new LTSA(data);
        expect(ltsa.parameter("neighbors")).toBe(2);
    });

    it("should throw error if d >= D", () => {
        const data = generateTestData(10, 2);
        expect(() => new LTSA(data, { d: 2 })).toThrow("Dimensionality of X (D = 2) must be greater than");
    });
});
