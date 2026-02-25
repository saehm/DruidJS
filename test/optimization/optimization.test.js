import { describe, expect, it } from "vitest";
import { powell } from "../../src/optimization/index.js";

describe("Optimization", () => {
    describe("powell", () => {
        it("should minimize a simple quadratic function", () => {
            const f = (x) => x[0] * x[0];
            const x0 = [5];
            const result = powell(f, x0);
            expect(result[0]).toBeCloseTo(0, 1);
        });

        it("should minimize a 2D quadratic function", () => {
            const f = (x) => x[0] * x[0] + x[1] * x[1];
            const x0 = [3, 4];
            const result = powell(f, x0);
            expect(result[0]).toBeCloseTo(0, 1);
            expect(result[1]).toBeCloseTo(0, 1);
        });

        it("should handle Rosenbrock function", () => {
            const f = (x) => (1 - x[0]) ** 2 + 100 * (x[1] - x[0] * x[0]) ** 2;
            const x0 = [-1, 2];
            const result = powell(f, x0, 1000);
            expect(result[0]).toBeCloseTo(1, 0);
            expect(result[1]).toBeCloseTo(1, 0);
        });

        it("should return Float64Array if input is Float64Array", () => {
            const f = (x) => x[0] * x[0];
            const x0 = new Float64Array([3]);
            const result = powell(f, x0);
            expect(result).toBeInstanceOf(Float64Array);
        });
    });
});
