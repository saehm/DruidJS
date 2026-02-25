import { describe, expect, it } from "vitest";
import { HierarchicalClustering } from "../../src/clustering/index.js";
import { euclidean } from "../../src/metrics/index.js";

describe("HierarchicalClustering", () => {
    it("should create hierarchical clusters with complete linkage", () => {
        const points = [
            [0, 0],
            [0, 1],
            [10, 10],
            [10, 11],
        ];
        const hc = new HierarchicalClustering(points, { linkage: "complete", metric: euclidean });

        expect(hc.root).toBeDefined();
        const clusters = hc.get_clusters(5, "distance");
        expect(clusters.length).toBeGreaterThan(0);

        const depthClusters = hc.get_clusters(1, "depth");
        expect(depthClusters.length).toBeGreaterThan(0);

        const rawDepthClusters = hc.get_clusters_raw(1, "depth");
        expect(rawDepthClusters.length).toBeGreaterThan(0);

        expect(() => hc.get_clusters(1, "invalid_type")).toThrow("invalid type");
    });

    it("should create hierarchical clusters with single linkage", () => {
        const points = [
            [0, 0],
            [0, 1],
            [10, 10],
            [10, 11],
        ];
        const hc = new HierarchicalClustering(points, { linkage: "single", metric: euclidean });

        const clusters = hc.get_clusters(5, "distance");
        expect(clusters.length).toBeGreaterThan(0);
    });

    it("should create hierarchical clusters with average linkage", () => {
        const points = [
            [0, 0],
            [0, 1],
            [10, 10],
            [10, 11],
        ];
        const hc = new HierarchicalClustering(points, { linkage: "average", metric: euclidean });

        const clusters = hc.get_clusters(5, "distance");
        expect(clusters.length).toBeGreaterThan(0);
    });

    it("should return cluster list", () => {
        const points = [
            [0, 0],
            [0, 1],
            [10, 10],
            [10, 11],
        ];
        const hc = new HierarchicalClustering(points, { linkage: "complete", metric: euclidean });

        const clusterList = hc.get_cluster_list(5, "distance");
        expect(clusterList).toHaveLength(4);
    });

    it("should handle invalid linkage gracefully", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        // Invalid linkage defaults to undefined behavior, but shouldn't crash
        const hc = new HierarchicalClustering(points, { linkage: "invalid", metric: euclidean });
        const clusters = hc.get_clusters(1, "distance");
        expect(clusters).toBeDefined();
    });

    it("should work with precomputed distance matrix", () => {
        // Provide a square distance matrix directly
        const distanceMatrix = [
            [0, 1, 10, 11],
            [1, 0, 11, 10],
            [10, 11, 0, 1],
            [11, 10, 1, 0],
        ];
        const hc = new HierarchicalClustering(distanceMatrix, {
            linkage: "complete",
            metric: "precomputed",
        });
        const clusters = hc.get_clusters(2, "distance");
        expect(clusters).toHaveLength(2);
    });

    it("should cut tree by depth", () => {
        const points = [
            [0, 0],
            [0, 1],
            [10, 10],
            [10, 11],
        ];
        const hc = new HierarchicalClustering(points, { linkage: "complete", metric: euclidean });
        const clusters = hc.get_clusters(1, "depth");
        expect(clusters.length).toBeGreaterThan(0);
    });

    it("should throw error for invalid cut type", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const hc = new HierarchicalClustering(points);
        // @ts-ignore
        expect(() => hc.get_clusters(1, "invalid")).toThrow("invalid type");
    });

    it("should return descendants for clusters", () => {
        const points = [
            [0, 0],
            [0, 1],
            [10, 10],
            [10, 11],
        ];
        const hc = new HierarchicalClustering(points);
        const rootDescendants = hc.root.descendants();
        expect(rootDescendants.length).toBeGreaterThan(0);
        expect(rootDescendants).toContain(hc.root);
    });
});
