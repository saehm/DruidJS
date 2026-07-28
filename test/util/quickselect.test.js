import { describe, expect, it } from "vitest";
import { quickselect, quickselectByAxis } from "../../src/util/quickselect.js";
import { Randomizer } from "../../src/util/randomizer.js";

/** Pivots come from a seeded randomizer, so every call below is reproducible. */
const R = new Randomizer(1212);

describe("quickselect", () => {
    it("should find the k-th smallest element in an unsorted array", () => {
        const arr = [9, 3, 7, 1, 5, 2, 8, 4, 6];
        // Sorted array: [1, 2, 3, 4, 5, 6, 7, 8, 9]

        // 0-th smallest -> 1
        expect(quickselect([...arr], R, 0)).toBe(1);
        // 4-th smallest (median) -> 5
        expect(quickselect([...arr], R, 4)).toBe(5);
        // 8-th smallest (max) -> 9
        expect(quickselect([...arr], R, 8)).toBe(9);
    });

    it("should work with custom comparison functions", () => {
        const items = [
            { id: "A", val: 50 },
            { id: "B", val: 10 },
            { id: "C", val: 30 },
            { id: "D", val: 20 },
            { id: "E", val: 40 },
        ];

        const medianItem = quickselect([...items], R, 2, (a, b) => a.val - b.val);
        expect(medianItem.val).toBe(30);
        expect(medianItem.id).toBe("C");
    });

    it("should handle edge cases (single element, empty array)", () => {
        expect(quickselect([42], R, 0)).toBe(42);
        expect(quickselect([], R, 0)).toBeUndefined();
    });
});

describe("quickselectByAxis", () => {
    it("should partition elements along a specified coordinate axis", () => {
        const points = [
            { index: 0, element: [10, 50] },
            { index: 1, element: [30, 10] },
            { index: 2, element: [20, 40] },
            { index: 3, element: [50, 20] },
            { index: 4, element: [40, 30] },
        ];

        // Partition along X-axis (axis 0)
        // X values: [10, 30, 20, 50, 40] -> sorted: 10, 20, 30, 40, 50
        const copyX = [...points];
        quickselectByAxis(copyX, R, 2, 0);
        expect(copyX[2].element[0]).toBe(30);

        // Partition along Y-axis (axis 1)
        // Y values: [50, 10, 40, 20, 30] -> sorted: 10, 20, 30, 40, 50
        const copyY = [...points];
        quickselectByAxis(copyY, R, 2, 1);
        expect(copyY[2].element[1]).toBe(30);
    });
});
