import { describe, expect, it } from "vitest";
import { Randomizer } from "../../src/util/randomizer.js";

describe("Randomizer", () => {
    it("should generate random numbers", () => {
        const r = new Randomizer(42);
        expect(r.random).toBeGreaterThanOrEqual(0);
        expect(r.random).toBeLessThanOrEqual(1);
        expect(r.seed).toBe(42);
    });

    it("should generate integers", () => {
        const r = new Randomizer(42);
        expect(Number.isInteger(r.random_int)).toBe(true);
    });

    it("should generate Gaussian random numbers", () => {
        const r = new Randomizer(42);
        const val1 = r.gauss_random();
        const val2 = r.gauss_random();
        expect(typeof val1).toBe("number");
        expect(typeof val2).toBe("number");
    });

    it("should handle choice from array", () => {
        const r = new Randomizer(42);
        const data = [1, 2, 3, 4, 5];
        const samples = r.choice(data, 2);
        expect(samples).toHaveLength(2);
        expect(data).toContain(samples[0]);
    });

    it("should throw error for invalid choice n", () => {
        const r = new Randomizer(42);
        expect(() => r.choice([1, 2], 5)).toThrow("n bigger than A!");
    });

    it("should throw error for non-array choice", () => {
        const r = new Randomizer(42);
        expect(() => r.choice("invalid", 2)).toThrow("A must be an Array!");
    });

    it("static choice works", () => {
        const samples = Randomizer.choice([1, 2, 3], 2, 123);
        expect(samples).toHaveLength(2);
    });
});
