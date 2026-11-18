import { expect } from "vitest";

/**
 * Validates that a result matrix/array contains only finite numbers
 * (no NaN, undefined, null, or Infinity).
 * @param {number[][] | import("../../src/matrix/index.js").Matrix} result
 * @param {string} [name="Result"] - Name of the algorithm for error messages
 */
export function expectValidValues(result, name = "Result") {
    const data = result.to2DArray ? result.to2DArray() : result;
    for (let i = 0; i < data.length; i++) {
        const row = data[i];
        for (let j = 0; j < row.length; j++) {
            const val = row[j];
            try {
                expect(Number.isFinite(val)).toBe(true);
                expect(val).not.toBeNaN();
                expect(val).not.toBeUndefined();
                expect(val).not.toBeNull();
            } catch (e) {
                console.error(`Validation failed for ${name} at [${i}, ${j}]: ${val}`);
                throw e;
            }
        }
    }
}

/**
 * The largest |Spearman rank correlation| of any column of `Y` with `target`. An embedding that
 * unrolls a manifold scores close to 1 against the manifold's own coordinates.
 *
 * @param {ArrayLike<ArrayLike<number>>} Y
 * @param {number[]} target
 * @returns {number}
 */
export function bestRankCorrelation(Y, target) {
    /** @param {number[]} c */
    const rank = (c) => {
        const r = new Array(c.length);
        c.map((v, i) => [v, i])
            .sort((a, b) => a[0] - b[0])
            .forEach(([, i], k) => {
                r[i] = k;
            });
        return r;
    };
    /**
     * @param {number[]} a
     * @param {number[]} b
     */
    const pearson = (a, b) => {
        const n = a.length;
        const ma = a.reduce((s, x) => s + x, 0) / n;
        const mb = b.reduce((s, x) => s + x, 0) / n;
        let ab = 0;
        let aa = 0;
        let bb = 0;
        for (let i = 0; i < n; ++i) {
            ab += (a[i] - ma) * (b[i] - mb);
            aa += (a[i] - ma) ** 2;
            bb += (b[i] - mb) ** 2;
        }
        return ab / Math.sqrt(aa * bb);
    };
    const rows = Array.from(Y, (row) => Array.from(row));
    const target_rank = rank(target);
    let best = 0;
    for (let k = 0; k < rows[0].length; ++k) {
        best = Math.max(best, Math.abs(pearson(rank(rows.map((r) => r[k])), target_rank)));
    }
    return best;
}
