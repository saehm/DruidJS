import { describe, expect, it } from "vitest";
import { minimum_spanning_tree } from "../../src/index.js";

/**
 * Total weight of an edge list.
 *
 * @param {[number, number, number][]} edges
 * @returns {number}
 */
const weight = (edges) => edges.reduce((sum, [, , w]) => sum + w, 0);

/**
 * Brute-force minimum spanning tree by trying every subset of `N - 1` edges, for checking the
 * result is genuinely optimal rather than merely a spanning tree.
 *
 * @param {[number, number, number][]} edges
 * @param {number} N
 * @returns {number}
 */
function brute_force_weight(edges, N) {
    let best = Infinity;
    const choose = (start, picked) => {
        if (picked.length === N - 1) {
            // connected?
            const parent = Array.from({ length: N }, (_, i) => i);
            const find = (x) => (parent[x] === x ? x : (parent[x] = find(parent[x])));
            let merges = 0;
            for (const [u, v] of picked) {
                const [a, b] = [find(u), find(v)];
                if (a !== b) {
                    parent[a] = b;
                    merges++;
                }
            }
            if (merges === N - 1) best = Math.min(best, weight(picked));
            return;
        }
        for (let i = start; i < edges.length; ++i) choose(i + 1, [...picked, edges[i]]);
    };
    choose(0, []);
    return best;
}

describe("minimum_spanning_tree", () => {
    /** @type {[number, number, number][]} */
    const diamond = [
        [0, 1, 1],
        [0, 2, 4],
        [1, 2, 2],
        [1, 3, 6],
        [2, 3, 3],
    ];

    it("should span every vertex with N-1 edges", () => {
        const tree = minimum_spanning_tree(diamond, 4);
        expect(tree).toHaveLength(3);
        const seen = new Set(tree.flatMap(([u, v]) => [u, v]));
        expect(seen.size).toBe(4);
    });

    it("should return the edges ascending by weight", () => {
        const tree = minimum_spanning_tree(diamond, 4);
        const weights = tree.map(([, , w]) => w);
        expect(weights).toEqual([...weights].sort((a, b) => a - b));
    });

    it("should find the actual minimum, not just some spanning tree", () => {
        expect(weight(minimum_spanning_tree(diamond, 4))).toBe(brute_force_weight(diamond, 4));
    });

    it("should find the minimum on a larger random graph", () => {
        let s = 99;
        const rand = () => {
            s = (s * 1103515245 + 12345) & 0x7fffffff;
            return s / 0x7fffffff;
        };
        const N = 7;
        /** @type {[number, number, number][]} */
        const edges = [];
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) edges.push([i, j, Math.round(rand() * 100)]);
        }
        expect(weight(minimum_spanning_tree(edges, N))).toBe(brute_force_weight(edges, N));
    });

    it("should not mutate the caller's edge list", () => {
        const input = diamond.map((e) => [...e]);
        const snapshot = JSON.stringify(input);
        minimum_spanning_tree(/** @type {any} */ (input), 4);
        expect(JSON.stringify(input)).toBe(snapshot);
    });

    it("should tolerate duplicate and reversed edges", () => {
        /** @type {[number, number, number][]} */
        const withDupes = [
            [0, 1, 1],
            [1, 0, 1],
            [0, 1, 5],
            [1, 2, 2],
            [2, 1, 9],
        ];
        const tree = minimum_spanning_tree(withDupes, 3);
        expect(tree).toHaveLength(2);
        expect(weight(tree)).toBe(3);
    });

    describe("degenerate input", () => {
        it("should return a forest when the graph is disconnected", () => {
            // two components: {0,1} and {2,3}, nothing bridging them
            /** @type {[number, number, number][]} */
            const split = [
                [0, 1, 1],
                [2, 3, 1],
            ];
            const tree = minimum_spanning_tree(split, 4);
            expect(tree).toHaveLength(2); // N-2, one edge short of spanning
        });

        it("should return nothing for fewer than two vertices", () => {
            expect(minimum_spanning_tree([], 0)).toEqual([]);
            expect(minimum_spanning_tree([], 1)).toEqual([]);
        });

        it("should return nothing when there are no edges", () => {
            expect(minimum_spanning_tree([], 5)).toEqual([]);
        });

        it("should reject an edge referring to a vertex outside [0, N)", () => {
            expect(() => minimum_spanning_tree([[0, 9, 1]], 3)).toThrow(/outside/);
            expect(() => minimum_spanning_tree([[-1, 0, 1]], 3)).toThrow(/outside/);
        });

        it("should include vertex 0 rather than treating it as absent", () => {
            // Regression: the disjoint set guarded with a falsy test, so vertex 0 -- which every
            // index-keyed caller starts at -- was reported as not found.
            const tree = minimum_spanning_tree(
                [
                    [0, 1, 1],
                    [1, 2, 1],
                ],
                3,
            );
            expect(tree).toHaveLength(2);
            expect(new Set(tree.flatMap(([u, v]) => [u, v])).has(0)).toBe(true);
        });
    });
});
