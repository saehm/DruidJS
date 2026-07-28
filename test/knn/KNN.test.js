import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { Annoy, BallTree, HNSW, KDTree, LSH, NaiveKNN, NNDescent } from "../../src/knn/index.js";
import { KNN } from "../../src/knn/KNN.js";
import { euclidean } from "../../src/metrics/index.js";

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

    it("should delegate search_by_index to the subclass search", () => {
        const knn = new KNN([[1, 2]]);
        // search_by_index is concrete now, so it fails through the still-abstract search.
        expect(() => knn.search_by_index(0, 5)).toThrow("The function search must be implemented!");
    });

    it("should return an empty result for an out-of-range index", () => {
        const knn = new KNN([[1, 2]]);
        expect(knn.search_by_index(-1, 5)).toEqual([]);
        expect(knn.search_by_index(99, 5)).toEqual([]);
    });
});

describe("KNN search_by_index contract", () => {
    // Every index must exclude the queried element itself. Callers rely on this instead of
    // stripping the self match themselves, which they each used to do differently.
    const points = mistle.IRIS.values.map((r) => Array.from(r));
    const k = 6;

    /** @type {[string, new (points: any, parameters?: any) => any][]} */
    const implementations = [
        ["Annoy", Annoy],
        ["BallTree", BallTree],
        ["HNSW", HNSW],
        ["KDTree", KDTree],
        ["LSH", LSH],
        ["NNDescent", NNDescent],
        ["NaiveKNN", NaiveKNN],
    ];

    for (const [name, Implementation] of implementations) {
        it(`${name} excludes the query element and returns k neighbors, closest first`, () => {
            const index = new Implementation(points, { metric: euclidean, seed: 42 });
            for (const i of [0, 17, 74, points.length - 1]) {
                const neighbors = index.search_by_index(i, k);
                expect(neighbors, `${name} returned nothing for ${i}`).toHaveLength(k);
                expect(neighbors.some((n) => n.index === i), `${name} returned self for ${i}`).toBe(false);
                expect(new Set(neighbors.map((n) => n.index)).size, `${name} returned a duplicate`).toBe(k);
                for (let j = 1; j < neighbors.length; ++j) {
                    expect(neighbors[j - 1].distance).toBeLessThanOrEqual(neighbors[j].distance);
                }
            }
        });
    }

    it("keeps duplicate points, excluding only the queried index", () => {
        // Three identical rows: querying one must still return the other two.
        const duplicated = [
            [1, 1],
            [1, 1],
            [1, 1],
            [9, 9],
        ];
        const index = new NaiveKNN(duplicated, { metric: euclidean, seed: 42 });
        const neighbors = index.search_by_index(0, 2);
        expect(neighbors.map((n) => n.index)).toEqual([1, 2]);
        expect(neighbors.every((n) => n.distance === 0)).toBe(true);
    });
});

describe("KNN parameter defaults", () => {
    // `KNN` stores the parameter object verbatim, so every subclass has to merge its own defaults.
    // Doing that with a default argument instead only covers the case where `parameters` is left
    // out entirely, so passing any single option silently dropped all the others. Every existing
    // test happened to pass the options it depended on, which left the default path unexercised —
    // NNDescent was throwing on construction the whole time.
    const points = mistle.IRIS.values;

    /** @type {[string, new (points: any, parameters?: any) => any][]} */
    const implementations = [
        ["Annoy", Annoy],
        ["BallTree", BallTree],
        ["HNSW", HNSW],
        ["KDTree", KDTree],
        ["LSH", LSH],
        ["NNDescent", NNDescent],
        ["NaiveKNN", NaiveKNN],
    ];

    for (const [name, Implementation] of implementations) {
        it(`${name} builds and searches with no parameters at all`, () => {
            const index = new Implementation(points);
            const neighbors = index.search(points[0], 5);

            expect(neighbors.length, `${name} returned nothing`).toBeGreaterThan(0);
            for (const { index: i, distance } of neighbors) {
                expect(i).toBeGreaterThanOrEqual(0);
                expect(i).toBeLessThan(points.length);
                expect(distance).toBeCloseTo(euclidean(points[i], points[0]), 12);
            }
        });

        it(`${name} keeps its remaining defaults when given one option`, () => {
            const index = new Implementation(points, { seed: 7 });

            expect(index._parameters.seed, `${name} ignored the option`).toBe(7);
            expect(index._parameters.metric, `${name} dropped its default metric`).toBe(euclidean);
            for (const [key, value] of Object.entries(index._parameters)) {
                expect(value, `${name} left '${key}' undefined`).toBeDefined();
            }

            expect(index.search(points[0], 5).length).toBeGreaterThan(0);
        });
    }
});
