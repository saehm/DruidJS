import { describe, expect, it } from "vitest";
import { KNN } from "../../src/knn/KNN.js";

describe("KNN Base Class", () => {
    it("should handle array elements", () => {
        const knn = new KNN([[1, 2]]);
        expect(knn._type).toBe("array");
    });

    it("should throw error for empty elements", () => {
        expect(() => new KNN([])).toThrow("Elements needs to contain at least one element!");
    });

    it("should handle typed array elements", () => {
        const knn = new KNN([new Float64Array([1, 2])]);
        expect(knn._type).toBe("typed");
    });

    it("should throw error on abstract method search", () => {
        const knn = new KNN([[1, 2]]);
        expect(() => knn.search([], 5)).toThrow("The function search must be implemented!");
    });

    it("should throw error on abstract method search_by_index", () => {
        const knn = new KNN([[1, 2]]);
        expect(() => knn.search_by_index(0, 5)).toThrow("The function search_by_index must be implemented!");
    });
});
