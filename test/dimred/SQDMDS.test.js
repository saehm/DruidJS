import { describe, expect, it } from "vitest";
import { SQDMDS } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("SQDMDS", () => {
    it("should complete within timeout", { timeout: 20000 }, () => {
        const data = generateTestData(20, 5);
        const sqdmds = new SQDMDS(data, { d: 2, seed: 42 });
        const result = sqdmds.transform(50);

        expect(result).toHaveLength(20);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 20000 }, () => {
        const data = generateTestData(15, 4);
        const sqdmds = new SQDMDS(data, { d: 2, seed: 42 });
        const result = sqdmds.transform(30);
        expectValidValues(result, "SQDMDS");
    });
    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = SQDMDS.transform(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = SQDMDS.generator(data, { d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await SQDMDS.transform_async(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should throw error if input is not square for precomputed", () => {
        const data = generateTestData(10, 5);
        expect(() => new SQDMDS(data, { metric: "precomputed" })).toThrow("SQDMDS input data must be a square Matrix");
    });
});
