import { describe, expect, it } from "vitest";
import { CURE } from "../../src/clustering/index.js";
import { euclidean, manhattan } from "../../src/metrics/index.js";

describe("CURE", () => {
    it("should cluster with default parameters", () => {
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
        const cure = new CURE(points, { K: 2, metric: euclidean });

        const clusters = cure.get_clusters();
        expect(clusters).toHaveLength(2);
    });

    it("should work with higher K", () => {
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
        const cure = new CURE(points, { K: 4, metric: euclidean });

        const clusters = cure.get_clusters();
        expect(clusters).toHaveLength(4);
    });

    it("should work with different metrics", () => {
        const points = [
            [0, 0],
            [3, 4],
            [6, 8],
        ];
        const cure = new CURE(points, { K: 2, metric: manhattan });

        const clusters = cure.get_clusters();
        expect(clusters).toHaveLength(2);
    });

    it("should work with custom num_representatives", () => {
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
        const cure = new CURE(points, { K: 2, num_representatives: 3, metric: euclidean });

        const clusters = cure.get_clusters();
        expect(clusters).toHaveLength(2);
    });

    it("should work with custom shrink_factor", () => {
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
        const cure = new CURE(points, { K: 2, shrink_factor: 0.3, metric: euclidean });

        const clusters = cure.get_clusters();
        expect(clusters).toHaveLength(2);
    });

    it("should handle K=1 (single cluster)", () => {
        const points = [
            [0, 0],
            [0, 1],
            [1, 0],
            [1, 1],
        ];
        const cure = new CURE(points, { K: 1, metric: euclidean });

        const clusters = cure.get_clusters();
        expect(clusters).toHaveLength(1);
        expect(clusters[0]).toHaveLength(4);
    });

    it("should handle K=N", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
        ];
        const cure = new CURE(points, { K: 3, metric: euclidean });

        const clusters = cure.get_clusters();
        expect(clusters).toHaveLength(3);
    });
});
