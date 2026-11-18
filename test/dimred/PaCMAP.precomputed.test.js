import { describe, expect, it } from "vitest";
import { LocalMAP, PaCMAP } from "../../src/dimred/index.js";
import { distance_matrix, Matrix } from "../../src/matrix/index.js";
import { euclidean, manhattan } from "../../src/metrics/index.js";
import { Randomizer } from "../../src/util/index.js";

/**
 * Fraction of each point's k nearest neighbors in `X` that stay among its k nearest in `Y`.
 *
 * @param {Matrix} X
 * @param {Matrix} Y
 * @param {number} k
 */
function neighbors_kept(X, Y, k = 10) {
    const DX = distance_matrix(X);
    const DY = distance_matrix(Y);
    const N = X.shape[0];
    /** @param {Matrix} D @param {number} i */
    const nn = (D, i) =>
        Array.from(D.row(i), (v, j) => [v, j])
            .filter(([, j]) => j !== i)
            .sort((a, b) => a[0] - b[0])
            .slice(0, k)
            .map(([, j]) => j);
    let kept = 0;
    for (let i = 0; i < N; ++i) {
        const near = new Set(nn(DX, i));
        for (const j of nn(DY, i)) if (near.has(j)) ++kept;
    }
    return kept / (N * k);
}

describe.each([
    ["PaCMAP", PaCMAP],
    ["LocalMAP", LocalMAP],
])('%s with metric "precomputed"', (_, Method) => {
    const r = new Randomizer(11);
    // three clusters in 6 dimensions; continuous values, so no ties between distances
    const X = new Matrix(150, 6, (i, j) => (j === i % 3 ? 4 : 0) + r.random);

    it.each([
        ["euclidean", euclidean],
        ["manhattan", manhattan],
    ])("builds the same pairs as from the data (%s)", { timeout: 60000 }, (_, metric) => {
        const direct = new Method(X, { metric, seed: 5 }).init();
        const pre = new Method(distance_matrix(X, metric), { metric: "precomputed", seed: 5 }).init();
        expect(pre._nn_pairs).toEqual(direct._nn_pairs);
        expect(pre._mn_pairs).toEqual(direct._mn_pairs);
        expect(pre._fp_pairs).toEqual(direct._fp_pairs);
    });

    it("keeps neighborhoods as well as from the data", { timeout: 60000 }, () => {
        const direct = /** @type {Matrix} */ (new Method(X, { seed: 5 }).transform());
        const pre = /** @type {Matrix} */ (
            new Method(distance_matrix(X), { metric: "precomputed", seed: 5 }).transform()
        );
        expect(pre.values.every(Number.isFinite)).toBe(true);
        expect(neighbors_kept(X, pre)).toBeGreaterThan(neighbors_kept(X, direct) - 0.03);
    });

    it("rejects a non-square matrix", () => {
        expect(() => new Method(new Matrix(20, 5, 1), { metric: "precomputed" }).init()).toThrow(/square/);
    });
});
