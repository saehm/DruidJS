import { describe, expect, it } from "vitest";
import { KMeans } from "../../src/clustering/index.js";
import { euclidean } from "../../src/metrics/index.js";

describe("KMeans", () => {
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
        const kmeans = new KMeans(points, { K: 2, metric: euclidean, seed: 42 });

        expect(kmeans.k).toBe(2);
        expect(kmeans.centroids).toHaveLength(2);

        const clusters = kmeans.get_clusters();
        expect(clusters).toHaveLength(2);
        expect(clusters[0].length + clusters[1].length).toBe(8);
    });

    it("should handle K larger than N by clamping K to N", () => {
        const points = [
            [0, 0],
            [10, 10],
        ];
        const kmeans = new KMeans(points, { K: 5, metric: euclidean });

        expect(kmeans.k).toBe(2);
    });

    it("should return cluster list", () => {
        const points = [
            [0, 0],
            [0, 1],
            [10, 10],
            [10, 11],
        ];
        const kmeans = new KMeans(points, { K: 2, metric: euclidean, seed: 42 });

        const clusterList = kmeans.get_cluster_list();
        expect(clusterList).toHaveLength(4);
        expect(new Set(clusterList).size).toBeLessThanOrEqual(2);
    });
});
