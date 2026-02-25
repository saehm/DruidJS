import { describe, expect, it } from "vitest";
import { kahan_sum, neumair_sum } from "../../src/numerical/index.js";

describe("Numerical Summation", () => {
    describe("kahan_sum", () => {
        it("should sum accurately", () => {
            const arr = [1, 2, 3, 4, 5];
            expect(kahan_sum(arr)).toBe(15);
        });

        it("should handle cancellation cases", () => {
            const arr = [1e16, 1, 1, 1, 1, -1e16];
            expect(kahan_sum(arr)).toBe(4);
        });
    });

    describe("neumair_sum", () => {
        it("should sum accurately", () => {
            const arr = [1, 2, 3, 4, 5];
            expect(neumair_sum(arr)).toBe(15);
        });

        it("should handle cancellation cases", () => {
            const arr = [1e16, 1, -1e16];
            expect(neumair_sum(arr)).toBe(1);
        });
    });
});
