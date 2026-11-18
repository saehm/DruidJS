import { describe, expect, it } from "vitest";
import { KMedoids } from "../../src/clustering/index.js";
import { distance_matrix, Matrix } from "../../src/matrix/index.js";
import { euclidean, manhattan } from "../../src/metrics/index.js";
import { Randomizer } from "../../src/util/index.js";

describe("KMedoids", () => {
    it("should cluster simple 2D data", () => {
        const points = [
            [0, 0],
            [0, 1],
            [1, 0],
            [1, 1],
            [10, 10],
            [10, 11],
            [11, 10],
            [11, 11],
        ];
        const kmedoids = new KMedoids(points, { K: 2, metric: euclidean, seed: 42 });

        const clusters = kmedoids.get_clusters();
        expect(clusters).toHaveLength(2);
        expect(clusters[0].length + clusters[1].length).toBe(8);
    });

    it("should return medoids", () => {
        const points = [
            [0, 0],
            [0, 1],
            [1, 0],
            [1, 1],
        ];
        const kmedoids = new KMedoids(points, { K: 2, metric: euclidean, seed: 42 });

        const medoids = kmedoids.get_medoids();
        expect(medoids).toBeDefined();
        expect(medoids.length).toBeGreaterThan(0);
    });

    it("should work with different metrics", () => {
        const points = [
            [0, 0],
            [3, 4],
            [6, 8],
        ];
        const kmedoids = new KMedoids(points, { K: 2, metric: manhattan });

        const clusters = kmedoids.get_clusters();
        expect(clusters).toHaveLength(2);
    });

    it("should handle K larger than N by clamping K to N", () => {
        const points = [
            [0, 0],
            [10, 10],
        ];
        const kmedoids = new KMedoids(points, { K: 5, metric: euclidean });

        expect(kmedoids.k).toBe(2);
    });
    it("should handle swaps during iteration", () => {
        // Points where random medoids are likely non-optimal
        const points = [
            [0, 0],
            [0, 1],
            [1, 0],
            [1, 1],
            [10, 10],
            [10, 11],
            [11, 10],
            [11, 11],
        ];
        const kmedoids = new KMedoids(points, { K: 2, seed: 123 });
        const clusters = kmedoids.get_clusters();
        expect(clusters.length).toBe(2);
    });

    it("should throw error if medoids are missing", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const kmedoids = new KMedoids(points, { K: 2 });
        // @ts-ignore
        kmedoids._cluster_medoids = [];
        // @ts-ignore
        expect(() => kmedoids._nearest_medoid([0, 0], 0)).toThrow("No medoids available");
    });

    it("should return the medoids and k", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const kmedoids = new KMedoids(points, { K: 2 });
        expect(kmedoids.medoids).toHaveLength(2);
        expect(kmedoids.k).toBe(2);
    });

    it.each([
        ["euclidean", euclidean],
        ["manhattan", manhattan],
    ])("gives the same clustering from a precomputed %s distance matrix", (_, metric) => {
        const r = new Randomizer(3);
        for (const N of [30, 200]) {
            const X = new Matrix(N, 3, () => r.random);
            const direct = new KMedoids(X, { K: 4, metric, seed: 7 });
            const pre = new KMedoids(distance_matrix(X, metric), { K: 4, metric: "precomputed", seed: 7 });
            expect(pre.get_medoids()).toEqual(direct.get_medoids());
            expect(pre.get_cluster_list()).toEqual(direct.get_cluster_list());
        }
    });

    it("rejects a non-square precomputed matrix", () => {
        expect(() => new KMedoids(new Matrix(4, 3, 1), { K: 2, metric: "precomputed" })).toThrow(/square/);
    });
});
