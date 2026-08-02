import { describe, expect, it } from "vitest";
import { MINFOTree } from "../../src/dimred/index.js";
import { generateClusteredData, generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

/**
 * Three well-separated gaussian blobs, so the tree has an unambiguous cluster structure to recover.
 *
 * @param {number} n_per - Points per blob.
 * @param {number} [seed]
 * @returns {{ X: number[][]; truth: number[] }}
 */
function blobs(n_per, seed = 42) {
    let s = seed;
    const rand = () => {
        s = (s * 1103515245 + 12345) & 0x7fffffff;
        return s / 0x7fffffff;
    };
    const gauss = () => Math.sqrt(-2 * Math.log(rand() + 1e-12)) * Math.cos(2 * Math.PI * rand());
    const centers = [
        [0, 0, 0],
        [6, 0, 0],
        [0, 6, 0],
    ];
    /** @type {number[][]} */
    const X = [];
    /** @type {number[]} */
    const truth = [];
    centers.forEach((center, ci) => {
        for (let i = 0; i < n_per; ++i) {
            X.push(center.map((c) => c + gauss() * 0.4));
            truth.push(ci);
        }
    });
    return { X, truth };
}

describe("MINFOTree", () => {
    it("should complete within timeout", { timeout: 10000 }, () => {
        const data = generateTestData(20, 4);
        const result = new MINFOTree(data, { clusters: 3, seed: 42 }).transform();

        expect(result).toHaveLength(20);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with generator", () => {
        const data = generateTestData(15, 5);
        const gen = new MINFOTree(data, { clusters: 2 }).generator();
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(15);
    });

    it("should work with static transform", () => {
        const data = generateTestData(15, 5);
        expect(MINFOTree.transform(data, { clusters: 2 })).toHaveLength(15);
    });

    it("should work with static generator", () => {
        const data = generateTestData(15, 5);
        const gen = MINFOTree.generator(data, { clusters: 2 });
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(15);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(15, 5);
        expect(await MINFOTree.transform_async(data, { clusters: 2 })).toHaveLength(15);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 4);
        const result = new MINFOTree(data, { clusters: 3, seed: 42 }).transform();
        expectValidValues(result, "MINFOTree");
    });

    it("should be deterministic for a given seed", () => {
        const data = generateTestData(20, 4);
        const a = new MINFOTree(data, { clusters: 3, seed: 42 }).transform();
        const b = new MINFOTree(data, { clusters: 3, seed: 42 }).transform();
        expect(a).toEqual(b);
    });

    it("should project to the requested dimensionality", () => {
        const data = generateTestData(20, 6);
        const result = new MINFOTree(data, { clusters: 3, d: 3, seed: 42 }).transform();
        expect(result[0]).toHaveLength(3);
    });

    describe("the spanning tree", () => {
        it("should be a spanning tree of all N points", () => {
            const { X } = blobs(15);
            const minfo = new MINFOTree(X, { clusters: 3, seed: 7 });
            const edges = minfo.edges;

            expect(edges).toHaveLength(X.length - 1);

            // every vertex is covered, and no edge is a self-loop
            const seen = new Set();
            for (const [u, v] of edges) {
                expect(u).not.toBe(v);
                seen.add(u);
                seen.add(v);
            }
            expect(seen.size).toBe(X.length);
        });

        it("should have strictly positive edge weights", () => {
            // A zero-weight edge would give the tree a zero-length path, and the 1/d² stiffness in
            // the Kamada-Kawai layout divides by it.
            const { X } = blobs(15);
            for (const [, , weight] of new MINFOTree(X, { clusters: 3, seed: 7 }).edges) {
                expect(weight).toBeGreaterThan(0);
            }
        });

        it("should keep clusters on separate branches", () => {
            // A tree over c clusters needs at least c-1 edges to join them. A faithful tree uses
            // close to that; a tree that scrambles the clusters uses many more.
            const { X, truth } = blobs(25);
            const minfo = new MINFOTree(X, { labels: truth, seed: 7 });
            const crossing = minfo.edges.filter(([u, v]) => truth[u] !== truth[v]).length;
            expect(crossing).toBeGreaterThanOrEqual(2);
            expect(crossing).toBeLessThanOrEqual(5);
        });

        it("should stay connected when the k-NN graph is not", () => {
            // Two far-apart clusters with k small enough that no k-NN edge bridges them.
            const data = generateClusteredData();
            const minfo = new MINFOTree(data, { clusters: 2, k: 2, seed: 42 });
            expect(minfo.edges).toHaveLength(data.length - 1);
            expectValidValues(minfo.transform(), "MINFOTree");
        });
    });

    describe("the Potts model", () => {
        it("should estimate a positive beta when labels agree with neighborhoods", () => {
            // Well-separated blobs: a point's neighbors almost always share its label, which is the
            // high-beta (ordered) regime of the Potts model.
            const { X, truth } = blobs(20);
            expect(new MINFOTree(X, { labels: truth, seed: 7 }).beta).toBeGreaterThan(1);
        });

        it("should estimate beta = 0 when labels are unrelated to the data", () => {
            // Random labels carry no spatial dependence, so observed energy cannot exceed the
            // expected energy at beta = 0 and the pseudo-likelihood equation has no positive root.
            const { X } = blobs(20);
            const labels = X.map((_, i) => i % 3);
            expect(new MINFOTree(X, { labels, seed: 7 }).beta).toBe(0);
        });

        it("should separate boundary points from interior points by curvature", () => {
            // Direction, not just separation: at the high beta well-separated clusters produce,
            // interior points end up at the TOP of the curvature range, not the bottom. That is the
            // opposite of the paper's prose in VI-C, and it follows from the formulas — phi and psi
            // both vanish for an interior point as beta grows, so `epsilon` in the denominator pins
            // S near 0. See MINFOTree._information_curvature for the measured crossover near beta 2.
            const { X, truth } = blobs(25);
            const minfo = new MINFOTree(X, { labels: truth, seed: 7 });
            const curvature = minfo.curvature;
            expect(minfo.beta).toBeGreaterThan(2);

            /** @type {Set<number>[]} */
            const neighborhood = Array.from({ length: X.length }, () => new Set());
            for (const [u, v] of minfo.edges) {
                neighborhood[u].add(v);
                neighborhood[v].add(u);
            }

            /** @type {number[]} */
            const boundary = [];
            /** @type {number[]} */
            const interior = [];
            for (let i = 0; i < X.length; ++i) {
                const mixes_labels = [...neighborhood[i]].some((j) => truth[j] !== truth[i]);
                (mixes_labels ? boundary : interior).push(curvature[i]);
            }

            expect(boundary.length).toBeGreaterThan(0);
            expect(interior.length).toBeGreaterThan(0);
            const mean = (/** @type {number[]} */ a) => a.reduce((x, y) => x + y, 0) / a.length;
            expect(mean(interior)).toBeGreaterThan(mean(boundary));
        });

        it("should recover the paper's curvature ordering at low beta", () => {
            // The regime is set by beta, not by epsilon. Corrupting a tenth of the labels loosens
            // the tie between a point and its neighborhood, which drops beta below 1 and restores
            // the ordering VI-C describes: interior points lower than boundary points.
            const { X, truth } = blobs(25);

            // deterministic 10% corruption
            let s = 99;
            const rand = () => {
                s = (s * 1103515245 + 12345) & 0x7fffffff;
                return s / 0x7fffffff;
            };
            const noisy = truth.map((t) => (rand() < 0.1 ? Math.floor(rand() * 3) : t));

            const minfo = new MINFOTree(X, { labels: noisy, seed: 7 });
            expect(minfo.beta).toBeLessThan(1);

            const curvature = minfo.curvature;
            /** @type {Set<number>[]} */
            const neighborhood = Array.from({ length: X.length }, () => new Set());
            for (const [u, v] of minfo.edges) {
                neighborhood[u].add(v);
                neighborhood[v].add(u);
            }
            /** @type {number[]} */
            const boundary = [];
            /** @type {number[]} */
            const interior = [];
            for (let i = 0; i < X.length; ++i) {
                const mixes_labels = [...neighborhood[i]].some((j) => noisy[j] !== noisy[i]);
                (mixes_labels ? boundary : interior).push(curvature[i]);
            }

            const mean = (/** @type {number[]} */ a) => a.reduce((x, y) => x + y, 0) / a.length;
            expect(mean(interior)).toBeLessThan(mean(boundary));
        });

        it("should normalize curvature into (0, 1]", () => {
            const { X } = blobs(15);
            for (const s of new MINFOTree(X, { clusters: 3, seed: 7 }).curvature) {
                expect(s).toBeGreaterThan(0);
                expect(s).toBeLessThanOrEqual(1 + 1e-3);
            }
        });
    });

    describe("the labels field", () => {
        it("should throw when neither clusters nor labels are given", () => {
            const data = generateTestData(15, 4);
            expect(() => new MINFOTree(data).transform()).toThrow(/labels field/);
        });

        it("should throw when labels length does not match X", () => {
            const data = generateTestData(15, 4);
            expect(() => new MINFOTree(data, { labels: [0, 1, 2] }).transform()).toThrow(/length/);
        });

        it("should accept non-numeric labels", () => {
            const { X, truth } = blobs(10);
            const named = truth.map((t) => ["alpha", "beta", "gamma"][t]);
            const minfo = new MINFOTree(X, { labels: named, seed: 7 });
            expect(new Set(minfo.labels).size).toBe(3);
            expectValidValues(minfo.transform(), "MINFOTree");
        });

        it("should recover well-separated clusters with Ward", () => {
            const { X, truth } = blobs(20);
            const found = new MINFOTree(X, { clusters: 3, seed: 7 }).labels;

            // Ward should split the blobs exactly, up to a relabeling.
            const mapping = new Map();
            for (let i = 0; i < X.length; ++i) {
                if (!mapping.has(truth[i])) mapping.set(truth[i], found[i]);
                expect(found[i]).toBe(mapping.get(truth[i]));
            }
            expect(new Set(mapping.values()).size).toBe(3);
        });

        it("should reject more clusters than a dendrogram cut can express", () => {
            // A dendrogram has N-1 merges and the traversal never emits a bare leaf, so N clusters
            // is unreachable. Cutting below every merge silently produces ONE cluster, so this has
            // to fail loudly rather than return the opposite of what was asked for.
            const { X } = blobs(6); // 18 points
            for (const clusters of [X.length, X.length + 5]) {
                expect(() => new MINFOTree(X, { clusters, seed: 7 }).transform()).toThrow(
                    /clusters must be at most/,
                );
            }
        });

        it("should accept the largest expressible cluster count", () => {
            const { X } = blobs(6);
            const minfo = new MINFOTree(X, { clusters: X.length - 1, seed: 7 });
            expect(new Set(minfo.labels).size).toBe(X.length - 1);
            expectValidValues(minfo.transform(), "MINFOTree");
        });

        it("should let kmeans go finer than the dendrogram can", () => {
            const { X } = blobs(6);
            const minfo = new MINFOTree(X, {
                clusters: X.length,
                clustering: "kmeans",
                seed: 7,
            });
            expect(new Set(minfo.labels).size).toBeGreaterThan(1);
            expectValidValues(minfo.transform(), "MINFOTree");
        });

        it("should support kmeans as the label source", () => {
            const { X } = blobs(15);
            const minfo = new MINFOTree(X, { clusters: 3, clustering: "kmeans", seed: 7 });
            expect(new Set(minfo.labels).size).toBe(3);
            expectValidValues(minfo.transform(), "MINFOTree");
        });
    });

    describe("the layout", () => {
        it("should reduce Kamada-Kawai energy below the MDS warm start", () => {
            const { X } = blobs(20);
            const minfo = new MINFOTree(X, { clusters: 3, seed: 7 });
            const edges = minfo.edges;
            const N = X.length;

            /** @type {[number, number][][]} */
            const adjacency = Array.from({ length: N }, () => []);
            for (const [u, v, w] of edges) {
                adjacency[u].push([v, w]);
                adjacency[v].push([u, w]);
            }
            const D = Array.from({ length: N }, () => new Float64Array(N));
            for (let source = 0; source < N; ++source) {
                const seen = new Uint8Array(N);
                seen[source] = 1;
                const stack = [source];
                while (stack.length) {
                    const u = /** @type {number} */ (stack.pop());
                    for (const [v, w] of adjacency[u]) {
                        if (seen[v]) continue;
                        seen[v] = 1;
                        D[source][v] = D[source][u] + w;
                        stack.push(v);
                    }
                }
            }

            /** @param {number[][]} Y */
            const energy = (Y) => {
                let E = 0;
                for (let i = 0; i < N; ++i) {
                    for (let j = i + 1; j < N; ++j) {
                        const residual = Math.hypot(Y[i][0] - Y[j][0], Y[i][1] - Y[j][1]) - D[i][j];
                        E += (residual * residual) / (D[i][j] * D[i][j]);
                    }
                }
                return E;
            };

            const mds = new MINFOTree(X, { clusters: 3, seed: 7, layout: "MDS" }).transform();
            const kk = new MINFOTree(X, { clusters: 3, seed: 7, layout: "kamada_kawai" }).transform();
            expect(energy(kk)).toBeLessThan(energy(mds));
        });

        it("should honour the MDS layout shortcut", () => {
            const { X } = blobs(15);
            const result = new MINFOTree(X, { clusters: 3, seed: 7, layout: "MDS" }).transform();
            expectValidValues(result, "MINFOTree");
        });
    });
});
