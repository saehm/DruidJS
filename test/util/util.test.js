import { describe, expect, it } from "vitest";
import { max } from "../../src/util/max.js";
import { min } from "../../src/util/min.js";
import { Randomizer } from "../../src/util/randomizer.js";

describe("Utils", () => {
    describe("max", () => {
        it("should return the maximum value", () => {
            expect(max([1, 5, 3, 9, 2])).toBe(9);
            expect(max([-1, -5, -2])).toBe(-1);
        });

        it("should handle null values", () => {
            expect(max([1, null, 5])).toBe(5);
        });
    });

    describe("min", () => {
        it("should return the minimum value", () => {
            expect(min([1, 5, 3, 9, 2])).toBe(1);
            expect(min([-1, -5, -2])).toBe(-5);
        });
    });

    describe("Randomizer", () => {
        it("should be deterministic with same seed", () => {
            const r1 = new Randomizer(42);
            const r2 = new Randomizer(42);
            for (let i = 0; i < 10; i++) {
                expect(r1.random).toBe(r2.random);
            }
        });

        it("should provide choice from array", () => {
            const r = new Randomizer(42);
            const arr = [1, 2, 3, 4, 5];
            const result = r.choice(arr, 3);
            expect(result).toHaveLength(3);
            result.forEach((item) => expect(arr).toContain(item));
        });
    });
});
