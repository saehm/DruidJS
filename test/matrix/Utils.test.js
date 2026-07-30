import { describe, expect, test } from "vitest";
import { distance_matrix, k_nearest_neighbors, linspace, Matrix, norm, normalize } from "../../src/matrix/index.js";
import { euclidean } from "../../src/metrics/index.js";

describe("Matrix Utilities", () => {
    describe("distance_matrix", () => {
        test("computes pairwise distances", () => {
            const data = Matrix.from([
                [0, 0],
                [3, 4],
            ]);
            const D = distance_matrix(data);
            expect(D.entry(0, 1)).toBe(5);
        });

        test("uses custom metric", () => {
            const data = Matrix.from([
                [0, 0],
                [3, 4],
            ]);
            const manhattan = (a, b) => Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]);
            const D = distance_matrix(data, manhattan);
            expect(D.entry(0, 1)).toBe(7);
        });
    });

    test("finds k nearest neighbors", () => {
        const points = [
            [0, 0],
            [1, 1],
            [10, 10],
        ];
        const knn = k_nearest_neighbors(points, 2, euclidean);
        expect(knn).toHaveLength(3);
        expect(knn[0]).toHaveLength(2);
        expect(knn[0][0].j).toBe(1);
        expect(knn[0][1].j).toBe(2);
    });

    describe("linspace", () => {
        test("generates evenly spaced values", () => {
            const result = linspace(0, 10, 5);
            expect(result).toEqual([0, 2.5, 5, 7.5, 10]);
        });

        test("defaults to one value per integer step when count is omitted", () => {
            expect(linspace(1, 5)).toEqual([1, 2, 3, 4, 5]);
            expect(linspace(0, 3)).toEqual([0, 1, 2, 3]);
        });

        test("treats an explicit null count the same as omitting it", () => {
            expect(linspace(1, 4, null)).toEqual(linspace(1, 4));
        });

        test("returns a single point when asked for one", () => {
            expect(linspace(2, 10, 1)).toEqual([2]);
        });

        test("returns empty when asked for none", () => {
            expect(linspace(0, 10, 0)).toEqual([]);
        });

        test("hits both endpoints exactly", () => {
            const result = linspace(-1, 1, 9);
            expect(result[0]).toBe(-1);
            expect(result.at(-1)).toBe(1);
            expect(result).toHaveLength(9);
        });

        test("counts down when end is before start", () => {
            expect(linspace(5, 1, 5)).toEqual([5, 4, 3, 2, 1]);
        });
    });

    describe("norm", () => {
        test("computes vector norm", () => {
            expect(norm([3, 4])).toBe(5);
        });
    });

    describe("normalize", () => {
        test("normalizes vector", () => {
            const v = [3, 4];
            const result = normalize(v);
            expect(norm(result)).toBeCloseTo(1, 10);
            expect(result[0]).toBeCloseTo(0.6, 10);
        });
    });
});
