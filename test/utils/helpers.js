import { expect } from "vitest";

/**
 * Validates that a result matrix/array contains only finite numbers
 * (no NaN, undefined, null, or Infinity).
 * @param {number[][] | import("../../src/matrix/index.js").Matrix} result
 * @param {string} [name="Result"] - Name of the algorithm for error messages
 */
export function expectValidValues(result, name = "Result") {
    const data = result.to2DArray ? result.to2DArray() : result;
    for (let i = 0; i < data.length; i++) {
        const row = data[i];
        for (let j = 0; j < row.length; j++) {
            const val = row[j];
            try {
                expect(Number.isFinite(val)).toBe(true);
                expect(val).not.toBeNaN();
                expect(val).not.toBeUndefined();
                expect(val).not.toBeNull();
            } catch (e) {
                console.error(`Validation failed for ${name} at [${i}, ${j}]: ${val}`);
                throw e;
            }
        }
    }
}
