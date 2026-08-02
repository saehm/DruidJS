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

    describe("ward linkage", () => {
        /** Three well-separated groups of four, so the correct 3-way split is unambiguous. */
        const blocks = [
            [0, 0],
            [0, 1],
            [1, 0],
            [1, 1],
            [20, 0],
            [20, 1],
            [21, 0],
            [21, 1],
            [0, 20],
            [0, 21],
            [1, 20],
            [1, 21],
        ];

        /** @param {HierarchicalClustering} hc */
        const merge_heights = (hc) =>
            hc.root
                .descendants()
                .filter((node) => !node.isLeaf)
                .map((node) => node.dist);

        it("should recover well-separated groups", () => {
            const hc = new HierarchicalClustering(blocks, { linkage: "ward", metric: euclidean });
            const heights = merge_heights(hc).sort((a, b) => a - b);
            // cut between the 3rd- and 2nd-tallest merge to leave exactly three clusters
            const cut = (heights[heights.length - 3] + heights[heights.length - 2]) / 2;
            const list = hc.get_cluster_list(cut);

            expect(new Set(list).size).toBe(3);
            // members of each block must agree, and differ from the other blocks
            for (const [a, b, c] of [
                [0, 1, 4],
                [4, 5, 8],
                [8, 9, 0],
            ]) {
                expect(list[a]).toBe(list[b]);
                expect(list[a]).not.toBe(list[c]);
            }
        });

        it("should produce monotone merge heights", () => {
            // Ward is a monotone (non-inverting) linkage. A mis-stated Lance-Williams recurrence
            // typically shows up as an inversion — a merge cheaper than one of its own children.
            const hc = new HierarchicalClustering(blocks, { linkage: "ward", metric: euclidean });
            /** @param {any} node */
            const check = (node) => {
                if (node.isLeaf) return;
                for (const child of [node.left, node.right]) {
                    if (!child || child.isLeaf) continue;
                    expect(node.dist).toBeGreaterThanOrEqual(child.dist - 1e-9);
                    check(child);
                }
            };
            check(hc.root);
        });

        it("should merge later than complete linkage on the same data", () => {
            // Ward's criterion grows with cluster size, so its dendrogram spans a wider range of
            // heights than a distance-based linkage on identical input.
            const ward = new HierarchicalClustering(blocks, { linkage: "ward", metric: euclidean });
            const complete = new HierarchicalClustering(blocks, {
                linkage: "complete",
                metric: euclidean,
            });
            expect(Math.max(...merge_heights(ward))).toBeGreaterThan(
                Math.max(...merge_heights(complete)),
            );
        });
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
