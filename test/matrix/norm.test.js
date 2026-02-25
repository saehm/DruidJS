import { describe, expect, it } from "vitest";
import { Matrix } from "../../src/matrix/index.js";
import { norm } from "../../src/matrix/norm.js";

describe("norm", () => {
    it("should compute norm from array", () => {
        expect(norm([1, 2, 3])).toBe(Math.sqrt(14));
    });

    it("should compute norm from a column matrix", () => {
        const v = Matrix.from([[1], [3]]);
        expect(norm(v)).toBe(Math.sqrt(10));
    });

    it("should compute norm from a row matrix", () => {
        const v = Matrix.from([[1, 3]]);
        expect(norm(v)).toBe(Math.sqrt(10));
    });

    it("should throw error for a 2d matrix", () => {
        const v = Matrix.from([
            [1, 2],
            [3, 4],
        ]);
        expect(() => norm(v)).toThrow("Matrix must be 1d!");
    });
});
