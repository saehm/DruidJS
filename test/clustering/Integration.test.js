import { describe, expect, it } from "vitest";
import {
    CURE,
    HierarchicalClustering,
    KMeans,
    KMedoids,
    MeanShift,
    OPTICS,
    XMeans,
} from "../../src/clustering/index.js";
import { euclidean } from "../../src/metrics/index.js";

describe("Clustering integration", () => {
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

    it("all clustering methods should work on same dataset", () => {
        const kmeans = new KMeans(points, { K: 2, metric: euclidean, seed: 42 });
        const kmedoids = new KMedoids(points, { K: 2, metric: euclidean, seed: 42 });
        const xmeans = new XMeans(points, { K_min: 2, K_max: 4, metric: euclidean, seed: 42 });
        const hc = new HierarchicalClustering(points, { linkage: "complete", metric: euclidean });
        const optics = new OPTICS(points, { epsilon: 5, min_points: 2, metric: euclidean });
        const meanShift = new MeanShift(points, { bandwidth: 3, metric: euclidean });
        const cure = new CURE(points, { K: 2, metric: euclidean });

        expect(kmeans.get_clusters().length).toBeGreaterThan(0);
        expect(kmedoids.get_clusters().length).toBeGreaterThan(0);
        expect(xmeans.get_clusters().length).toBeGreaterThan(0);
        expect(hc.get_clusters(10, "distance").length).toBeGreaterThan(0);
        expect(optics.get_clusters().length).toBeGreaterThan(0);
        expect(meanShift.get_clusters().length).toBeGreaterThan(0);
        expect(cure.get_clusters().length).toBeGreaterThan(0);
    });

    it("clustering results should partition all points", () => {
        const kmeans = new KMeans(points, { K: 2, metric: euclidean, seed: 42 });
        const clusters = kmeans.get_clusters();
        const totalPoints = clusters.reduce((sum, cluster) => sum + cluster.length, 0);
        expect(totalPoints).toBe(8);
    });
});
