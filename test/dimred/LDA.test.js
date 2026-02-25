import { describe, expect, it } from "vitest";
import { LDA } from "../../src/dimred/index.js";
import { generateClusteredData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("LDA", () => {
    it("should reduce dimensionality to specified d", () => {
        const data = generateClusteredData();
        const labels = Array(5).fill(0).concat(Array(5).fill(1));
        const lda = new LDA(data, { labels, d: 1, seed: 42 });
        const result = lda.transform();

        expect(result).toHaveLength(10);
        expect(result[0]).toHaveLength(1);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateClusteredData();
        const labels = Array(5).fill(0).concat(Array(5).fill(1));
        const lda = new LDA(data, { labels, d: 1, seed: 42 });
        const result = lda.transform();
        expectValidValues(result, "LDA");
    });
    it("should work with static transform", () => {
        const data = generateClusteredData();
        const labels = Array(5).fill(0).concat(Array(5).fill(1));
        const result = LDA.transform(data, { labels, d: 1 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateClusteredData();
        const labels = Array(5).fill(0).concat(Array(5).fill(1));
        const gen = LDA.generator(data, { labels, d: 1 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateClusteredData();
        const labels = Array(5).fill(0).concat(Array(5).fill(1));
        const result = await LDA.transform_async(data, { labels, d: 1 });
        expect(result).toHaveLength(10);
    });

    it("should throw error if labels are missing or length mismatch", () => {
        const data = generateClusteredData();
        const lda = new LDA(data, { labels: null });
        expect(() => lda.transform()).toThrow("LDA needs parameter label");

        const lda2 = new LDA(data, { labels: [0] });
        expect(() => lda2.transform()).toThrow("LDA needs parameter label");
    });

    it("should use random seed in eig_args if not provided", () => {
        const data = generateClusteredData();
        const labels = Array(5).fill(0).concat(Array(5).fill(1));
        const lda = new LDA(data, { labels });
        // @ts-ignore
        expect(lda.parameter("eig_args").seed).toBeDefined();
    });
});
