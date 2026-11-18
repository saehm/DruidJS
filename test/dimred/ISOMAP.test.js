import * as mistle from "@saehrimnir/mistle";
import { describe, expect, it } from "vitest";
import { ISOMAP, MDS } from "../../src/dimred/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("ISOMAP", () => {
    it("should use the largest positive eigenvalues on a disconnected graph", () => {
        // A ring of 8 points and a separate triangle: with 2 neighbors the graph has two components,
        // and Isomap fills in the missing distances with 10 × the largest geodesic distance. Geodesics
        // along the ring are not Euclidean, so the double-centered matrix has the eigenvalues
        // 2039.84, 8, 8, 1.734, ..., -2.343, -2.343 (numpy). By magnitude, -2.343 would take the
        // fourth axis, which then had no extent; it belongs to 1.734.
        /** @type {number[][]} */
        const X = [];
        for (let i = 0; i < 8; ++i) {
            const t = (2 * Math.PI * i) / 8 + 0.1;
            X.push([Math.cos(t), Math.sin(t)]);
        }
        X.push([50, 0], [51, 0], [50, 1.5]);

        const Y = /** @type {number[][]} */ (new ISOMAP(X, { d: 4, neighbors: 2 }).transform());
        expectValidValues(Y, "ISOMAP");
        const spread = (/** @type {number} */ k) => {
            const c = Y.map((row) => row[k]);
            const m = c.reduce((s, x) => s + x, 0) / c.length;
            return Math.sqrt(c.reduce((s, x) => s + (x - m) ** 2, 0) / c.length);
        };
        // sqrt(λ / n) for the four largest eigenvalues by value
        const expected = [13.6176249424, 0.8528028654, 0.8528028654, 0.3970507714];
        for (let k = 0; k < 4; ++k) expect(spread(k)).toBeCloseTo(expected[k], 5);
    });

    it("should equal classical MDS when every point is a neighbor", () => {
        // With the complete kNN graph, geodesic distances are the Euclidean distances, so the
        // embedding (including the sqrt(eigenvalue) scaling of each axis) must match MDS.
        const X = mistle.IRIS.values;
        const I = /** @type {number[][]} */ (new ISOMAP(X, { d: 2, neighbors: X.length - 1 }).transform());
        const M = /** @type {number[][]} */ (new MDS(X, { d: 2 }).transform());
        for (let i = 0; i < X.length; ++i) {
            for (let k = 0; k < 2; ++k) {
                expect(Math.abs(I[i][k])).toBeCloseTo(Math.abs(M[i][k]), 6);
            }
        }
    });

    it("should reduce dimensionality to d=2 by default", () => {
        const data = generateTestData(15, 4);
        const isomap = new ISOMAP(data, { neighbors: 4, seed: 42 });
        const result = isomap.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 4);
        const isomap = new ISOMAP(data, { neighbors: 4, d: 2, seed: 42 });
        const result = isomap.transform();
        expectValidValues(result, "ISOMAP");
    });
    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = ISOMAP.transform(data, { neighbors: 4, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 4);
        const gen = ISOMAP.generator(data, { neighbors: 2, d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 4);
        const result = await ISOMAP.transform_async(data, { neighbors: 2, d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should use default neighbors if not provided", () => {
        const data = generateTestData(20, 5);
        const isomap = new ISOMAP(data);
        expect(isomap.parameter("neighbors")).toBe(2);
    });

    it("should work with SMACOF projection", () => {
        const data = generateTestData(15, 4);
        const isomap = new ISOMAP(data, { neighbors: 4, seed: 42, project: "SMACOF" });
        const result = isomap.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });
});
