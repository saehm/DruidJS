import { afterEach, describe, expect, it } from "vitest";
import { KKMDS, Matrix, SAMMON, setWasmEnabled, SMACOF } from "../../src/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

/**
 * A swiss roll, which has enough curvature that the weighting scheme visibly matters.
 *
 * @param {number} n
 * @param {number} [seed]
 * @returns {number[][]}
 */
function swissroll(n, seed = 7) {
    let s = seed;
    const rand = () => {
        s = (s * 1103515245 + 12345) & 0x7fffffff;
        return s / 0x7fffffff;
    };
    /** @type {number[][]} */
    const X = [];
    for (let i = 0; i < n; ++i) {
        const t = 1.5 * Math.PI * (1 + 2 * rand());
        X.push([t * Math.cos(t), 21 * rand(), t * Math.sin(t)]);
    }
    return X;
}

/**
 * Fraction of each point's `k` true nearest neighbors that survive into the embedding.
 *
 * @param {number[][]} X
 * @param {number[][]} Y
 * @param {number} [k]
 * @returns {number}
 */
function neighborhood_preservation(X, Y, k = 10) {
    /**
     * @param {number[][]} P
     * @param {number} i
     * @returns {number[]}
     */
    const nearest = (P, i) =>
        P.map((p, j) => [j, Math.hypot(...p.map((v, c) => v - P[i][c]))])
            .filter(([j]) => j !== i)
            .sort((a, b) => a[1] - b[1])
            .slice(0, k)
            .map(([j]) => j);

    let hits = 0;
    for (let i = 0; i < X.length; ++i) {
        const truth = new Set(nearest(X, i));
        for (const j of nearest(Y, i)) if (truth.has(j)) hits++;
    }
    return hits / (X.length * k);
}

describe("KKMDS", () => {
    afterEach(() => setWasmEnabled(true));

    it("should complete within timeout", { timeout: 10000 }, () => {
        const data = generateTestData(20, 4);
        const result = new KKMDS(data, { seed: 42 }).transform();
        expect(result).toHaveLength(20);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with generator", () => {
        const data = generateTestData(20, 5);
        const gen = new KKMDS(data, { seed: 42 }).generator();
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(20);
    });

    it("should work with static transform", () => {
        expect(KKMDS.transform(generateTestData(20, 5), { seed: 42 })).toHaveLength(20);
    });

    it("should work with static generator", () => {
        const gen = KKMDS.generator(generateTestData(20, 5), { seed: 42 });
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(20);
    });

    it("should work with static transform_async", async () => {
        expect(await KKMDS.transform_async(generateTestData(20, 5), { seed: 42 })).toHaveLength(20);
    });

    it("should produce valid values", () => {
        expectValidValues(new KKMDS(generateTestData(20, 4), { seed: 42 }).transform(), "KKMDS");
    });

    it("should be deterministic for a given seed", () => {
        const data = generateTestData(25, 4);
        expect(new KKMDS(data, { seed: 42 }).transform()).toEqual(new KKMDS(data, { seed: 42 }).transform());
    });

    it("should project to the requested dimensionality", () => {
        const result = new KKMDS(generateTestData(20, 6), { d: 3, seed: 42 }).transform();
        expect(result[0]).toHaveLength(3);
    });

    it("should return a Matrix for Matrix input", () => {
        const Y = new KKMDS(Matrix.from(generateTestData(20, 4)), { seed: 42 }).transform();
        expect(Y).toBeInstanceOf(Matrix);
        expect(Y.shape).toEqual([20, 2]);
    });

    describe("the objective", () => {
        it("should decrease the energy monotonically", () => {
            const X = swissroll(40);
            const kkmds = new KKMDS(X, { seed: 42 });
            kkmds.check_init();
            const initial = kkmds.energy;

            // The line search only ever accepts a strictly lower energy, so no step may raise it.
            let previous = Infinity;
            let steps = 0;
            for (const _ of kkmds.generator()) {
                expect(kkmds.energy).toBeLessThanOrEqual(previous);
                previous = kkmds.energy;
                steps++;
            }
            expect(steps).toBeGreaterThan(1);
            expect(kkmds.energy).toBeLessThan(initial);
        });

        it("should always yield at least the warm start", () => {
            // Even when the line search rejects the first step, the generator must produce a
            // configuration — otherwise `transform` would have nothing to return.
            const X = generateTestData(20, 4);
            const yielded = [...new KKMDS(X, { seed: 42, iterations: 0 }).generator()];
            expect(yielded.length).toBeGreaterThanOrEqual(1);
            expect(yielded[0]).toHaveLength(20);
        });

        it("should beat a random start from the MDS warm start", () => {
            // The whole reason init_DR defaults to "MDS": the energy is non-convex, and Kamada-Kawai
            // is classically sensitive to where the descent begins.
            const X = swissroll(60);
            const warm = new KKMDS(X, { seed: 42, init_DR: "MDS" });
            warm.transform();
            const cold = new KKMDS(X, { seed: 42, init_DR: "random" });
            cold.transform();
            expect(warm.energy).toBeLessThan(cold.energy);
        });

        it("should preserve neighborhoods better than SMACOF and SAMMON as configured", () => {
            // Deliberately an as-configured comparison of the three classes, not an attribution to
            // the 1/d^2 weighting. Holding the solver fixed and varying only the exponent produces a
            // much smaller and data-dependent effect -- see the exponent tests in StressMDS.test.js.
            // What this pins is the user-facing claim: reaching for KKMDS over SMACOF or SAMMON on
            // manifold data does preserve more of the neighborhood structure.
            const X = swissroll(120);
            const kkmds = neighborhood_preservation(X, new KKMDS(X, { seed: 42 }).transform());
            const smacof = neighborhood_preservation(X, new SMACOF(X, { seed: 42 }).transform());
            const sammon = neighborhood_preservation(
                X,
                new SAMMON(X, { seed: 42, init_DR: "MDS" }).transform(),
            );

            expect(kkmds).toBeGreaterThan(smacof);
            expect(kkmds).toBeGreaterThan(sammon);
        });
    });

    describe("robustness", () => {
        it("should accept a precomputed distance matrix", () => {
            const X = swissroll(30);
            const D = new Matrix(X.length, X.length, (i, j) =>
                Math.hypot(...X[i].map((v, c) => v - X[j][c])),
            );
            expectValidValues(new KKMDS(D, { metric: "precomputed", seed: 42 }).transform(), "KKMDS");
        });

        it("should tolerate duplicate rows", () => {
            // A duplicate gives d_ij = 0 and so an infinite 1/d^2 stiffness. Those pairs are dropped.
            const X = swissroll(25);
            X.push([...X[0]], [...X[0]], [...X[3]]);
            expectValidValues(new KKMDS(X, { seed: 42 }).transform(), "KKMDS");
        });

        it("should handle zero iterations as the warm start alone", () => {
            const X = swissroll(20);
            const none = new KKMDS(X, { seed: 42, iterations: 0 }).transform();
            const mds_only = new KKMDS(X, { seed: 42, init_DR: "MDS", iterations: 0 }).transform();
            expect(none).toEqual(mds_only);
        });

        it("should reject PCA init on a precomputed matrix", () => {
            const D = new Matrix(10, 10, (i, j) => Math.abs(i - j));
            expect(() => new KKMDS(D, { metric: "precomputed", init_DR: "PCA" }).transform()).toThrow(
                /precomputed/,
            );
        });

        it("should reject an unknown init_DR", () => {
            expect(() =>
                new KKMDS(generateTestData(20, 4), { init_DR: /** @type {any} */ ("nope") }).transform(),
            ).toThrow(/init_DR/);
        });
    });

    describe("the WASM kernel", () => {
        it("should agree with the JS implementation", () => {
            // Above WASM_MIN_ROWS the kernel takes over; the two paths must not diverge.
            const X = swissroll(80);

            setWasmEnabled(true);
            const accelerated = new KKMDS(X, { seed: 42 });
            const Y_wasm = accelerated.transform();

            setWasmEnabled(false);
            const plain = new KKMDS(X, { seed: 42 });
            const Y_js = plain.transform();

            expect(Y_wasm).toHaveLength(Y_js.length);
            for (let i = 0; i < Y_js.length; ++i) {
                for (let k = 0; k < Y_js[i].length; ++k) {
                    expect(Y_wasm[i][k]).toBeCloseTo(Y_js[i][k], 6);
                }
            }
            expect(accelerated.energy).toBeCloseTo(plain.energy, 4);
        });

        it("should release its buffers when a run ends", async () => {
            const { held_session_count } = await import("../../src/wasm/index.js");
            setWasmEnabled(true);
            new KKMDS(swissroll(40), { seed: 42 }).transform();
            expect(held_session_count()).toBe(0);
        });

        it("should release its buffers when a generator is abandoned", async () => {
            const { held_session_count } = await import("../../src/wasm/index.js");
            setWasmEnabled(true);
            const kkmds = new KKMDS(swissroll(40), { seed: 42 });
            const steps = kkmds.generator();
            steps.next();
            kkmds.release();
            expect(held_session_count()).toBe(0);
        });
    });
});
