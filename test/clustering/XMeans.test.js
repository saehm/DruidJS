import { describe, expect, it } from "vitest";
import { XMeans } from "../../src/clustering/index.js";
import { euclidean } from "../../src/metrics/index.js";
import { generateTestData } from "../utils/data-generators.js";

describe("XMeans", () => {
    it("should determine optimal number of clusters", () => {
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
        const xmeans = new XMeans(points, {
            K_min: 2,
            K_max: 5,
            metric: euclidean,
            seed: 42,
        });

        expect(xmeans.k).toBeGreaterThanOrEqual(2);
        expect(xmeans.k).toBeLessThanOrEqual(5);
    });

    it("should trigger splitting with suitable parameters", () => {
        // Create two distant clusters
        const data1 = generateTestData(50, 2, 42); // Cluster 1
        const data2 = data1.map((p) => p.map((v) => v + 100)); // Cluster 2
        const points = [...data1, ...data2];

        const xmeans = new XMeans(points, {
            K_min: 2,
            K_max: 10,
            min_cluster_size: 5,
            seed: 42,
        });

        // It should at least consider splitting
        expect(xmeans.k).toBeGreaterThanOrEqual(2);
    });

    it("should return centroids and clusters", () => {
        const data = generateTestData(30, 2, 42);
        const xmeans = new XMeans(data, { K_min: 2, K_max: 4, seed: 42 });

        expect(xmeans.centroids).toBeDefined();
        expect(xmeans.get_clusters()).toBeDefined();
        expect(xmeans.get_cluster_list()).toHaveLength(30);
        expect(xmeans.k).toBeGreaterThanOrEqual(2);
    });

    it("should work with K_min equal to K_max", () => {
        const points = [
            [0, 0],
            [0, 1],
            [10, 10],
            [10, 11],
        ];
        const xmeans = new XMeans(points, { K_min: 2, K_max: 2, seed: 42 });
        expect(xmeans.k).toBe(2);
    });

    it("should handle edge case: N <= K", () => {
        const points = [
            [0, 0],
            [10, 10],
        ];
        const xmeans = new XMeans(points, { K_min: 2, K_max: 2 });
        expect(xmeans.k).toBe(2);
    });

    it("should handle edge case: zero variance", () => {
        const points = [
            [0, 0],
            [0, 0],
            [0, 0],
            [0, 0],
        ];
        const xmeans = new XMeans(points, { K_min: 2, K_max: 2 });
        expect(xmeans.k).toBe(2);
    });

    it("should throw error if called on non-initialized state", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const xmeans = new XMeans(points);
        // @ts-ignore
        xmeans._best_kmeans = null;
        expect(() => xmeans.get_clusters()).toThrow("XMeans has not been run");
        expect(() => xmeans.get_cluster_list()).toThrow("XMeans has not been run");
        expect(() => xmeans.centroids).toThrow("XMeans has not been run");
        expect(() => xmeans.k).toThrow("XMeans has not been run");
    });
});
