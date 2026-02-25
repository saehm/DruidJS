import { describe, expect, it } from "vitest";
import { SMACOF } from "../../src/dimred/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("SMACOF", () => {
    it("should reduce dimensionality to specified d", () => {
        const data = generateTestData(15, 4);
        const smacof = new SMACOF(data, { d: 2, seed: 42 });
        const result = smacof.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 5);
        const smacof = new SMACOF(data, { d: 2, seed: 42 });
        const result = smacof.transform();
        expectValidValues(result, "SMACOF");
    });

    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = SMACOF.transform(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = SMACOF.generator(data, { d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await SMACOF.transform_async(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with precomputed metric", () => {
        const data = generateTestData(10, 5);
        const dists = Matrix.from(generateTestData(10, 10)); // dummy dist matrix
        const smacof = new SMACOF(dists, { metric: "precomputed", d: 2 });
        const result = smacof.transform();
        expect(result.shape[0]).toBe(10);
    });
});
