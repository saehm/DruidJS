import { describe, expect, it } from "vitest";
import { MeanShift } from "../../src/clustering/index.js";
import { euclidean, manhattan } from "../../src/metrics/index.js";

describe("MeanShift", () => {
    it("should cluster with flat kernel", () => {
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
        const meanShift = new MeanShift(points, { kernel: "flat", bandwidth: 3, metric: euclidean });

        const clusters = meanShift.get_clusters();
        expect(clusters.length).toBeGreaterThan(0);
    });

    it("should work with gaussian kernel", () => {
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
        const meanShift = new MeanShift(points, {
            kernel: "gaussian",
            bandwidth: 3,
            metric: euclidean,
        });

        const clusters = meanShift.get_clusters();
        expect(clusters.length).toBeGreaterThan(0);
    });

    it("should work with custom kernel function", () => {
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
        const customKernel = (dist) => Math.exp(-dist);
        const meanShift = new MeanShift(points, {
            kernel: customKernel,
            bandwidth: 3,
            metric: euclidean,
        });

        const clusters = meanShift.get_clusters();
        expect(clusters.length).toBeGreaterThan(0);
    });

    it("should return cluster list", () => {
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
        const meanShift = new MeanShift(points, { metric: euclidean });

        const clusterList = meanShift.get_cluster_list();
        expect(clusterList).toHaveLength(8);
        for (const id of clusterList) {
            expect(id).toBeGreaterThanOrEqual(0);
        }
    });

    it("should work with different metrics", () => {
        const points = [
            [0, 0],
            [3, 4],
            [6, 8],
        ];
        const meanShift = new MeanShift(points, { bandwidth: 5, metric: manhattan });

        const clusters = meanShift.get_clusters();
        expect(clusters.length).toBeGreaterThan(0);
    });

    it("should handle lazy initialization if clusters are cleared", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const meanShift = new MeanShift(points);
        // @ts-ignore
        meanShift._cluster_list = undefined;
        // @ts-ignore
        meanShift._clusters = undefined;

        const clusters = meanShift.get_clusters();
        expect(clusters.length).toBeGreaterThan(0);

        const list = meanShift.get_cluster_list();
        expect(list.length).toBe(2);
    });
});
