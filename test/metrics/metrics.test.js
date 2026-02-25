import { describe, expect, it } from "vitest";
import {
    bray_curtis,
    canberra,
    chebyshev,
    cosine,
    euclidean,
    euclidean_squared,
    goodman_kruskal,
    hamming,
    haversine,
    jaccard,
    manhattan,
    sokal_michener,
    wasserstein,
    yule,
} from "../../src/metrics/index.js";

describe("euclidean", () => {
    it("should compute euclidean distance between two vectors", () => {
        const a = [0, 0];
        const b = [3, 4];
        expect(euclidean(a, b)).toBe(5);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 2, 3];
        const b = [1, 2, 3];
        expect(euclidean(a, b)).toBe(0);
    });

    it("should work with Float64Array", () => {
        const a = new Float64Array([0, 0]);
        const b = new Float64Array([3, 4]);
        expect(euclidean(a, b)).toBe(5);
    });

    it("should throw error for different length vectors", () => {
        expect(() => euclidean([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });

    it("should handle negative values", () => {
        const a = [-1, -2];
        const b = [2, 3];
        expect(euclidean(a, b)).toBe(Math.sqrt(9 + 25));
    });
});

describe("euclidean_squared", () => {
    it("should compute squared euclidean distance", () => {
        const a = [0, 0];
        const b = [3, 4];
        expect(euclidean_squared(a, b)).toBe(25);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 2, 3];
        const b = [1, 2, 3];
        expect(euclidean_squared(a, b)).toBe(0);
    });

    it("should throw error for different length vectors", () => {
        expect(() => euclidean_squared([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });
});

describe("manhattan", () => {
    it("should compute manhattan distance", () => {
        const a = [0, 0];
        const b = [3, 4];
        expect(manhattan(a, b)).toBe(7);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 2, 3];
        const b = [1, 2, 3];
        expect(manhattan(a, b)).toBe(0);
    });

    it("should handle negative values", () => {
        const a = [-1, -2];
        const b = [2, 3];
        expect(manhattan(a, b)).toBe(8);
    });

    it("should throw error for different length vectors", () => {
        expect(() => manhattan([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });
});

describe("chebyshev", () => {
    it("should compute chebyshev distance", () => {
        const a = [0, 0];
        const b = [3, 4];
        expect(chebyshev(a, b)).toBe(4);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 2, 3];
        const b = [1, 2, 3];
        expect(chebyshev(a, b)).toBe(0);
    });

    it("should handle negative values", () => {
        const a = [-5, 2];
        const b = [3, -8];
        expect(chebyshev(a, b)).toBe(10);
    });

    it("should throw error for different length vectors", () => {
        expect(() => chebyshev([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });
});

describe("cosine", () => {
    it("should compute cosine distance (angle in radians)", () => {
        const a = [1, 0];
        const b = [1, 1];
        expect(cosine(a, b)).toBeCloseTo(Math.PI / 4, 10);
    });

    it("should return 0 for parallel vectors", () => {
        const a = [1, 2];
        const b = [2, 4];
        expect(cosine(a, b)).toBeCloseTo(0, 5);
    });

    it("should return π for opposite vectors", () => {
        const a = [1, 0];
        const b = [-1, 0];
        expect(cosine(a, b)).toBeCloseTo(Math.PI, 10);
    });

    it("should throw error for different length vectors", () => {
        expect(() => cosine([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });
});

describe("jaccard", () => {
    it("should compute jaccard distance for binary-like vectors", () => {
        const a = [1, 0, 1, 0, 1];
        const b = [1, 0, 0, 1, 1];
        expect(jaccard(a, b)).toBe(0.5);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 0, 1];
        const b = [1, 0, 1];
        expect(jaccard(a, b)).toBe(0);
    });

    it("should return 1 for completely disjoint vectors", () => {
        const a = [1, 0, 0];
        const b = [0, 1, 1];
        expect(jaccard(a, b)).toBe(1);
    });

    it("should throw error for different length vectors", () => {
        expect(() => jaccard([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });
});

describe("hamming", () => {
    it("should compute hamming distance", () => {
        const a = [1, 0, 1, 0];
        const b = [1, 1, 0, 0];
        expect(hamming(a, b)).toBe(0.5);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 2, 3];
        const b = [1, 2, 3];
        expect(hamming(a, b)).toBe(0);
    });

    it("should return 1 for completely different vectors", () => {
        const a = [1, 2];
        const b = [3, 4];
        expect(hamming(a, b)).toBe(1);
    });

    it("should throw error for different length vectors", () => {
        expect(() => hamming([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });
});

describe("canberra", () => {
    it("should compute canberra distance", () => {
        const a = [0, 1, 2];
        const b = [1, 0, 3];
        expect(canberra(a, b)).toBe(2.2);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 2, 3];
        const b = [1, 2, 3];
        expect(canberra(a, b)).toBe(0);
    });

    it("should handle zeros by returning NaN when both are zero", () => {
        const a = [0, 1];
        const b = [0, 2];
        expect(canberra(a, b)).toBeNaN();
    });

    it("should throw error for different length vectors", () => {
        expect(() => canberra([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });
});

describe("bray_curtis", () => {
    it("should compute bray-curtis distance", () => {
        const a = [1, 2, 3];
        const b = [4, 5, 6];
        expect(bray_curtis(a, b)).toBeCloseTo(9 / 21, 10);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 2, 3];
        const b = [1, 2, 3];
        expect(bray_curtis(a, b)).toBe(0);
    });

    it("should return 1 for completely disjoint vectors", () => {
        const a = [1, 0, 0];
        const b = [0, 1, 1];
        expect(bray_curtis(a, b)).toBe(1);
    });

    it("should throw error for different length vectors", () => {
        expect(() => bray_curtis([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });
});

describe("haversine", () => {
    it("should compute haversine distance on unit sphere", () => {
        const a = [0, 0];
        const b = [0, Math.PI / 2];
        expect(haversine(a, b)).toBeCloseTo(Math.PI / 2, 10);
    });

    it("should return 0 for identical points", () => {
        const a = [0.5, 0.5];
        const b = [0.5, 0.5];
        expect(haversine(a, b)).toBe(0);
    });

    it("should throw error for wrong dimensions", () => {
        expect(() => haversine([0, 0, 0], [0, 0])).toThrow("Haversine distance requires exactly 2 coordinates");
        expect(() => haversine([0, 0], [0, 0, 0])).toThrow("Haversine distance requires exactly 2 coordinates");
    });
});

describe("wasserstein", () => {
    it("should compute wasserstein (EMD) distance", () => {
        const a = [1, 0, 0];
        const b = [0, 0, 1];
        expect(wasserstein(a, b)).toBe(2);
    });

    it("should return 0 for identical distributions", () => {
        const a = [0.25, 0.25, 0.25, 0.25];
        const b = [0.25, 0.25, 0.25, 0.25];
        expect(wasserstein(a, b)).toBe(0);
    });

    it("should throw error for different length vectors", () => {
        expect(() => wasserstein([1, 2], [1, 2, 3])).toThrow("Vector a and b needs to be of the same length!");
    });

    it("should return 0 if both distributions sum to 0", () => {
        expect(wasserstein([0, 0], [0, 0])).toBe(0);
    });

    it("should return Infinity if only one distribution sums to 0", () => {
        expect(wasserstein([0, 0], [1, 1])).toBe(Infinity);
        expect(wasserstein([1, 1], [0, 0])).toBe(Infinity);
    });
});

describe("yule", () => {
    it("should compute yule distance", () => {
        const a = [1, 0, 1, 0];
        const b = [1, 1, 0, 0];
        expect(yule(a, b)).toBe(1);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 0, 1];
        const b = [1, 0, 1];
        expect(yule(a, b)).toBe(0);
    });
});

describe("sokal_michener", () => {
    it("should compute sokal-michener distance", () => {
        const a = [1, 0, 1, 0];
        const b = [1, 1, 0, 0];
        expect(sokal_michener(a, b)).toBe(2 / 3);
    });

    it("should return 0 for identical vectors", () => {
        const a = [1, 0, 1];
        const b = [1, 0, 1];
        expect(sokal_michener(a, b)).toBe(0);
    });
});

describe("goodman_kruskal", () => {
    it("should compute goodman-kruskal gamma coefficient", () => {
        const a = [1, 2, 3, 4];
        const b = [1, 2, 3, 4];
        expect(goodman_kruskal(a, b)).toBe(1);
    });

    it("should return -1 for perfect inverse correlation", () => {
        const a = [1, 2, 3, 4];
        const b = [4, 3, 2, 1];
        expect(goodman_kruskal(a, b)).toBe(-1);
    });

    it("should handle ties in only one sequence", () => {
        const a = [1, 1, 2, 3]; // tie in a
        const b = [1, 2, 3, 4]; // no tie in b
        expect(goodman_kruskal(a, b)).not.toBeNaN();

        const c = [1, 2, 3, 4]; // no tie in c
        const d = [1, 1, 2, 3]; // tie in d
        expect(goodman_kruskal(c, d)).not.toBeNaN();
    });
});
