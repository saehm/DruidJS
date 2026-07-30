import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { LSH } from "../../src/knn/index.js";
import { euclidean } from "../../src/metrics/index.js";

describe("LSH", () => {
    it("should find nearest neighbors in 2D", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
            [10, 10],
        ];
        const lsh = new LSH(points, { metric: euclidean });

        const neighbors = lsh.search([0, 0], 2);
        expect(neighbors).toHaveLength(2);
        expect(neighbors[0].index).toBe(0);
    });

    it("should work with IRIS dataset", () => {
        const points = mistle.IRIS.values;
        const lsh = new LSH(points, { metric: euclidean });

        const neighbors = lsh.search(points[0], 3);
        expect(neighbors).toHaveLength(3);
        expect(neighbors[0].index).toBe(0);
    });

    it("should handle large dataset (YEAST)", () => {
        const points = mistle.YEAST.values;
        const lsh = new LSH(points, { metric: euclidean });
        const neighbors = lsh.search(points[0], 10);

        expect(neighbors).toHaveLength(10);
    });

    it("should work with search_by_index", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
        ];
        const lsh = new LSH(points, { metric: euclidean });
        const neighbors = lsh.search_by_index(0, 2);
        expect(neighbors).toHaveLength(2);
        // search_by_index excludes the queried element itself.
        expect(neighbors.map((n) => n.index)).toEqual([1, 2]);
    });

    it("should fallback to linear search when candidates are insufficient", () => {
        const points = [
            [0, 0],
            [1, 1],
            [2, 2],
        ];
        // Use 1 hash table to increase chance of insufficient candidates for large k
        const lsh = new LSH(points, { metric: euclidean, numHashTables: 1 });
        const neighbors = lsh.search([0, 0], 5);
        expect(neighbors.length).toBeLessThanOrEqual(3);
    });

    describe("built empty and filled with add()", () => {
        const points = mistle.YEAST.values;

        /**
         * Counts how many distances a run of queries actually computes. This is what separates a
         * working index from one that answers correctly by scanning everything — recall alone
         * cannot tell them apart, because brute force scores a perfect 1.0.
         *
         * @param {(metric: (a: any, b: any) => number) => LSH<any>} build
         */
        function distances_per_query(build) {
            let calls = 0;
            const counting = (/** @type {any} */ a, /** @type {any} */ b) => {
                calls++;
                return euclidean(a, b);
            };
            const lsh = build(counting);
            const queries = points.slice(0, 20);
            calls = 0;
            for (const q of queries) lsh.search(q, 10);
            return calls / queries.length;
        }

        it("should build a real index rather than one bucket per table", () => {
            // `_dim` and the bucket width both come from the data. Derived from the placeholder
            // `super` gets, the projection vectors are length 1, `_computeHash` reads past their
            // end, and every point hashes to "NaN,NaN,…" — one bucket per table, every query a
            // full scan.
            const lsh = new LSH([], { metric: euclidean, seed: 1212 });
            lsh.add(points);

            const buckets = lsh._hashTables.reduce((sum, table) => sum + table.size, 0);
            expect(buckets).toBeGreaterThan(lsh._hashTables.length);
            for (const table of lsh._hashTables) {
                for (const hash of table.keys()) expect(hash).not.toContain("NaN");
            }
        });

        it("should examine a small fraction of the dataset per query", () => {
            const incremental = distances_per_query((metric) => {
                const lsh = new LSH([], { metric, seed: 1212 });
                lsh.add(points);
                return lsh;
            });
            expect(incremental).toBeLessThan(points.length / 2);
        });

        it("should match an index constructed with the same points directly", () => {
            const direct = new LSH(points, { metric: euclidean, seed: 1212 });
            const incremental = new LSH([], { metric: euclidean, seed: 1212 });
            incremental.add(points);

            expect(incremental._dim).toBe(direct._dim);
            expect(incremental._bucketWidth).toBe(direct._bucketWidth);
            for (const q of points.slice(0, 10)) {
                expect(incremental.search(q, 5).map((n) => n.index)).toEqual(direct.search(q, 5).map((n) => n.index));
            }
        });

        it("should find the query itself when it is in the index", () => {
            const lsh = new LSH([], { metric: euclidean, seed: 1212 });
            lsh.add(points);
            for (const i of [0, 17, 128]) {
                expect(lsh.search(points[i], 5)[0].index).toBe(i);
            }
        });

        it("should keep indices contiguous across several add() calls", () => {
            const lsh = new LSH([], { metric: euclidean, seed: 1212 });
            lsh.add(points.slice(0, 100));
            lsh.add(points.slice(100));
            expect(lsh._elements).toHaveLength(points.length);
            expect(lsh.search(points[250], 1)[0].index).toBe(250);
        });

        it("should respect an explicitly supplied bucket width", () => {
            const lsh = new LSH([], { metric: euclidean, seed: 1212, bucketWidth: 7 });
            lsh.add(points);
            expect(lsh._bucketWidth).toBe(7);
        });

        it("should search an index that was never filled", () => {
            expect(new LSH([], { metric: euclidean, seed: 1212 }).search([1, 2, 3], 3)).toEqual([]);
        });
    });
});
