import { describe, expect, it } from "vitest";
import { OPTICS } from "../../src/clustering/index.js";
import { euclidean } from "../../src/metrics/index.js";

describe("OPTICS", () => {
    it("should cluster with density-based algorithm", () => {
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
        const optics = new OPTICS(points, { epsilon: 2, min_points: 2, metric: euclidean });

        const clusters = optics.get_clusters();
        expect(clusters.length).toBeGreaterThan(0);
    });

    it("should return cluster list", () => {
        const points = [
            [0, 0],
            [0, 1],
            [1, 0],
            [1, 1],
            [10, 10],
        ];
        const optics = new OPTICS(points, { epsilon: 2, min_points: 2, metric: euclidean });

        const clusterList = optics.get_cluster_list();
        expect(clusterList).toHaveLength(5);
    });

    it("should identify outliers", () => {
        const points = [
            [0, 0],
            [0.1, 0.1],
            [0.2, 0.2],
            [100, 100],
        ];
        const optics = new OPTICS(points, { epsilon: 0.5, min_points: 2, metric: euclidean });

        const clusterList = optics.get_cluster_list();
        expect(clusterList).toContain(-1);
    });

    it("should work with different epsilon values", () => {
        const points = [
            [0, 0],
            [0, 1],
            [1, 0],
            [1, 1],
        ];
        const optics = new OPTICS(points, { epsilon: 1, min_points: 2, metric: euclidean });

        const clusters = optics.get_clusters();
        expect(clusters.length).toBeGreaterThan(0);
    });

    it("should handle single cluster case", () => {
        const points = [
            [0, 0],
            [0.1, 0.1],
            [0.2, 0.2],
        ];
        const optics = new OPTICS(points, { epsilon: 1, min_points: 2, metric: euclidean });

        const clusters = optics.get_clusters();
        expect(clusters.length).toBeGreaterThan(0);
    });
});
