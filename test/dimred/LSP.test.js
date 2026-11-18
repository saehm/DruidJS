import { describe, expect, it } from "vitest";
import { LSP } from "../../src/dimred/index.js";
import { distance_matrix, Matrix } from "../../src/matrix/index.js";
import { euclidean, manhattan } from "../../src/metrics/index.js";
import { generateTestData } from "../utils/data-generators.js";

describe("LSP", () => {
    it("should reduce dimensionality to specified d", () => {
        const data = generateTestData(30, 5);
        const lsp = new LSP(data, { neighbors: 5, control_points: 10, d: 2, seed: 42 });
        const result = lsp.transform();

        expect(result).toHaveLength(30);
        expect(result[0]).toHaveLength(2);
    });

    it("should build a Laplacian whose rows sum to zero", () => {
        // Each row is 1 on the diagonal plus exactly `neighbors` entries of -1/neighbors. The
        // neighbor search used to include the point itself and then skip it, which left only
        // neighbors - 1 entries and a row sum of 1/neighbors.
        const K = 5;
        const N = 30;
        const data = generateTestData(N, 5);
        const lsp = new LSP(data, { neighbors: K, control_points: 10, d: 2, seed: 42 });
        lsp.transform();
        // @ts-ignore - private
        const A = lsp._A;

        for (let i = 0; i < N; ++i) {
            let sum = 0;
            let off_diagonal = 0;
            for (let j = 0; j < N; ++j) {
                sum += A.entry(i, j);
                if (i !== j && A.entry(i, j) !== 0) ++off_diagonal;
            }
            expect(off_diagonal, `row ${i}`).toBe(K);
            expect(sum, `row ${i}`).toBeCloseTo(0, 12);
        }
    });

    it("should use default parameters if none provided", () => {
        const data = generateTestData(20, 5);
        const lsp = new LSP(data);
        expect(lsp.parameter("neighbors")).toBeGreaterThan(0);
        expect(lsp.parameter("control_points")).toBeGreaterThan(0);
    });

    it("should handle double initialization", () => {
        const data = generateTestData(10, 5);
        const lsp = new LSP(data);
        lsp.init();
        const firstA = lsp._A;
        lsp.init();
        expect(lsp._A).toBe(firstA);
    });

    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = LSP.transform(data, { d: 2 });
        expect(result).toHaveLength(10);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = LSP.generator(data, { d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await LSP.transform_async(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it.each([
        ["euclidean", euclidean],
        ["manhattan", manhattan],
    ])("gives the same projection from a precomputed %s distance matrix", (_, metric) => {
        const X = Matrix.from(generateTestData(60, 5));
        const direct = /** @type {Matrix} */ (new LSP(X, { metric, seed: 3 }).transform());
        const pre = /** @type {Matrix} */ (
            new LSP(distance_matrix(X, metric), { metric: "precomputed", seed: 3 }).transform()
        );
        let diff = 0;
        for (let i = 0; i < direct.values.length; ++i)
            diff = Math.max(diff, Math.abs(direct.values[i] - pre.values[i]));
        expect(diff).toBeLessThan(1e-8);
    });

    it("rejects a non-square precomputed matrix", () => {
        const lsp = new LSP(new Matrix(10, 3, 1), { metric: "precomputed" });
        expect(() => lsp.transform()).toThrow(/square/);
    });
});
