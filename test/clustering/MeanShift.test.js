import { afterEach, beforeEach, describe, expect, it } from "vitest";
import { MeanShift } from "../../src/clustering/index.js";
import { euclidean, manhattan } from "../../src/metrics/index.js";
import { setWasmEnabled } from "../../src/wasm/index.js";

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

    describe("synchronous update (JS path)", () => {
        /**
         * Mean shift moves every point using the positions held at the *start* of the iteration.
         * `Matrix.row` returns a live subarray, so writing the shift straight back into the matrix
         * would let point `i` see the already-moved `0..i-1` — turning the result into a function
         * of row order and, worse, disagreeing with the WASM kernel, which reads all of its
         * neighbours from the untouched input.
         *
         * WASM is switched off for the whole block: with the kernel available these never reach the
         * JS loop, and would pass no matter what it does.
         */
        beforeEach(() => setWasmEnabled(false));
        afterEach(() => setWasmEnabled(true));

        const cloud = [
            [0, 0],
            [0.2, 0.1],
            [0.1, 0.25],
            [5, 5],
            [5.2, 5.1],
            [5.1, 5.25],
            [2.4, 2.6],
        ];

        /**
         * Where each point ended up, keyed by the coordinates it started from — so the result can
         * be compared across input orderings. Converged positions are the sharp test: the final
         * partition can survive an order-dependent sweep, the modes generally do not.
         *
         * @param {number[][]} points
         * @param {object} opts
         */
        function final_positions(points, opts) {
            const ms = new MeanShift(points, { seed: 42, ...opts });
            ms.get_cluster_list();
            return new Map(points.map((p, i) => [JSON.stringify(p), Array.from(ms._points.row(i))]));
        }

        /**
         * @param {number[][]} points
         * @param {object} opts
         */
        function order_sensitivity(points, opts) {
            const forward = final_positions(points, opts);
            const reversed = final_positions([...points].reverse(), opts);
            let worst = 0;
            for (const [key, p] of forward) {
                const q = /** @type {number[]} */ (reversed.get(key));
                for (let d = 0; d < p.length; ++d) worst = Math.max(worst, Math.abs(p[d] - q[d]));
            }
            return worst;
        }

        for (const opts of [
            { bandwidth: 1, kernel: "flat" },
            { bandwidth: 1, kernel: "gaussian" },
            { bandwidth: 2, kernel: "gaussian" },
        ]) {
            it(`should converge to the same modes whatever order the points arrive in — ${opts.kernel}, bandwidth ${opts.bandwidth}`, () => {
                // An in-place sweep moves these by 1e-1 or so when the rows are reversed.
                expect(order_sensitivity(cloud, opts)).toBeLessThan(1e-12);
            });
        }

        it("should place a point independently of the points listed before it", () => {
            // Same cloud, one point relocated to the front. Every *other* point must land exactly
            // where it did before — it cannot, if earlier rows' shifts are visible to later ones.
            const moved_to_front = [cloud[6], ...cloud.slice(0, 6)];
            const before = final_positions(cloud, { bandwidth: 1, kernel: "gaussian" });
            const after = final_positions(moved_to_front, { bandwidth: 1, kernel: "gaussian" });

            for (const [key, p] of before) {
                const q = /** @type {number[]} */ (after.get(key));
                for (let d = 0; d < p.length; ++d) expect(q[d], `${key}[${d}]`).toBeCloseTo(p[d], 12);
            }
        });
    });

    it("should handle lazy initialization if clusters are cleared", () => {
        const points = [
            [0, 0],
            [1, 1],
        ];
        const meanShift = new MeanShift(points);
        // @ts-expect-error
        meanShift._cluster_list = undefined;
        // @ts-expect-error
        meanShift._clusters = undefined;

        const clusters = meanShift.get_clusters();
        expect(clusters.length).toBeGreaterThan(0);

        const list = meanShift.get_cluster_list();
        expect(list.length).toBe(2);
    });
});
