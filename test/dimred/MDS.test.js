import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { MDS, PCA } from "../../src/dimred/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

/**
 * @param {number[][]} Y
 * @param {number} k
 */
const column = (Y, k) => Y.map((row) => row[k]);

/** @param {number[]} c */
const std = (c) => {
    const m = c.reduce((s, x) => s + x, 0) / c.length;
    return Math.sqrt(c.reduce((s, x) => s + (x - m) ** 2, 0) / c.length);
};

/**
 * @param {number[]} a
 * @param {number[]} b
 */
const correlation = (a, b) => {
    const ma = a.reduce((s, x) => s + x, 0) / a.length;
    const mb = b.reduce((s, x) => s + x, 0) / b.length;
    let ab = 0;
    for (let i = 0; i < a.length; ++i) ab += (a[i] - ma) * (b[i] - mb);
    return ab / a.length / (std(a) * std(b));
};

describe("MDS", () => {
    // Classical MDS on Euclidean distances is PCA up to the sign of each axis. Iris and Wine have
    // distinct leading eigenvalues, so the axes are unique.
    it.each([
        ["Iris", mistle.IRIS.values],
        ["Wine", mistle.WINE.values],
    ])("should equal PCA up to sign on %s", (_, X) => {
        const P = /** @type {number[][]} */ (new PCA(X, { d: 2 }).transform());
        const M = /** @type {number[][]} */ (new MDS(X, { d: 2 }).transform());
        for (let k = 0; k < 2; ++k) {
            expect(Math.abs(correlation(column(M, k), column(P, k)))).toBeCloseTo(1, 6);
            expect(std(column(M, k)) / std(column(P, k))).toBeCloseTo(1, 6);
        }
    });

    it("should give no extent to axes with negative eigenvalues", () => {
        // Not Euclidean (d(2, 3) = 3 > d(2, 0) + d(0, 3) = 2): the double-centered matrix has the
        // eigenvalues 4.5, -1.5, 0.5 and 0. Each axis has spread sqrt(max(λ, 0) / n).
        const D = [
            [0, 1, 1, 1],
            [1, 0, 1, 1],
            [1, 1, 0, 3],
            [1, 1, 3, 0],
        ];
        const Y = /** @type {Matrix} */ (new MDS(Matrix.from(D), { d: 3, metric: "precomputed" }).transform());
        const extents = [0, 1, 2].map((k) => std(Array.from(Y.col(k))));
        expect(extents[0]).toBeCloseTo(Math.sqrt(4.5 / 4), 6);
        expect(extents[1]).toBeCloseTo(0, 6);
        expect(extents[2]).toBeCloseTo(Math.sqrt(0.5 / 4), 6);
        expectValidValues(Y.to2dArray(), "MDS");
    });

    it("should reduce dimensionality to specified d", () => {
        const data = generateTestData(15, 4);
        const mds = new MDS(data, { d: 2, seed: 42 });
        const result = mds.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 5);
        const mds = new MDS(data, { d: 2, seed: 42 });
        const result = mds.transform();
        expectValidValues(result, "MDS");
    });
    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = MDS.transform(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = MDS.generator(data, { d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await MDS.transform_async(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with precomputed metric", () => {
        const data = generateTestData(10, 5);
        const dists = Matrix.from(generateTestData(10, 10)); // dummy dist matrix
        const mds = new MDS(dists, { metric: "precomputed", d: 2 });
        const result = mds.transform();
        expect(result.shape[0]).toBe(10);
    });

    it("should compute stress", () => {
        const data = generateTestData(10, 5);
        const mds = new MDS(data, { d: 2 });
        mds.transform();
        const stress = mds.stress();
        expect(stress).toBeGreaterThanOrEqual(0);
    });

    it("should throw error if stress is called before transform", () => {
        const data = generateTestData(10, 5);
        const mds = new MDS(data, { d: 2 });
        expect(() => mds.stress()).toThrow("First transform!");
    });
});
