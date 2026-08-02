import { afterEach, describe, expect, it } from "vitest";
import {
    KKMDS,
    Matrix,
    SAMMON,
    setWasmEnabled,
    StressMDS,
    WEIGHTS_ELASTIC,
    WEIGHTS_SAMMON,
    WEIGHTS_UNIFORM,
} from "../../src/index.js";
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
 * The objective itself, computed independently of the implementation under test.
 *
 * @param {number[][]} Y
 * @param {number[][]} X
 * @param {number} q
 * @returns {number}
 */
function weighted_stress(Y, X, q) {
    let total = 0;
    for (let i = 0; i < X.length; ++i) {
        for (let j = i + 1; j < X.length; ++j) {
            const d = Math.hypot(...X[i].map((v, c) => v - X[j][c]));
            if (!(d > 0)) continue;
            const r = Math.hypot(...Y[i].map((v, c) => v - Y[j][c])) - d;
            total += d ** q * r * r;
        }
    }
    return total;
}

describe("StressMDS", () => {
    afterEach(() => setWasmEnabled(true));

    it("should complete within timeout", { timeout: 10000 }, () => {
        const result = new StressMDS(generateTestData(20, 4), { seed: 42 }).transform();
        expect(result).toHaveLength(20);
        expect(result[0]).toHaveLength(2);
    });

    it("should work with generator", () => {
        const gen = new StressMDS(generateTestData(20, 5), { seed: 42 }).generator();
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(20);
    });

    it("should work with static transform", () => {
        expect(StressMDS.transform(generateTestData(20, 5), { seed: 42 })).toHaveLength(20);
    });

    it("should work with static generator", () => {
        const gen = StressMDS.generator(generateTestData(20, 5), { seed: 42 });
        let result;
        for (const val of gen) {
            result = val;
        }
        expect(result).toHaveLength(20);
    });

    it("should work with static transform_async", async () => {
        expect(await StressMDS.transform_async(generateTestData(20, 5), { seed: 42 })).toHaveLength(20);
    });

    it("should produce valid values", () => {
        expectValidValues(new StressMDS(generateTestData(20, 4), { seed: 42 }).transform(), "StressMDS");
    });

    it("should be deterministic for a given seed", () => {
        const data = generateTestData(25, 4);
        expect(new StressMDS(data, { seed: 42 }).transform()).toEqual(
            new StressMDS(data, { seed: 42 }).transform(),
        );
    });

    it("should return a Matrix for Matrix input", () => {
        const Y = new StressMDS(Matrix.from(generateTestData(20, 4)), { seed: 42 }).transform();
        expect(Y).toBeInstanceOf(Matrix);
        expect(Y.shape).toEqual([20, 2]);
    });

    describe("the weighting", () => {
        it("should minimize the objective it was given", () => {
            // The real contract: an exponent q must produce the lowest sigma_q of any exponent. This
            // is checked against a stress function written independently of the implementation.
            const X = swissroll(60);
            const solutions = new Map(
                [WEIGHTS_UNIFORM, WEIGHTS_SAMMON, WEIGHTS_ELASTIC].map((q) => [
                    q,
                    new StressMDS(X, { weights: q, seed: 42 }).transform(),
                ]),
            );

            for (const q of solutions.keys()) {
                const own = weighted_stress(/** @type {number[][]} */ (solutions.get(q)), X, q);
                for (const other of solutions.keys()) {
                    if (other === q) continue;
                    const rival = weighted_stress(/** @type {number[][]} */ (solutions.get(other)), X, q);
                    expect(own).toBeLessThan(rival);
                }
            }
        });

        it("should trade global fidelity away as the exponent drops", () => {
            // The one effect of the exponent that is monotone and solver-independent: concentrating
            // the objective on short distances costs global fidelity. Verified monotone in 9/9
            // size/seed combinations, so the ordering is asserted in full rather than at the ends.
            //
            // The converse -- that a lower exponent buys *local* fidelity -- is deliberately NOT
            // asserted: it holds only on data with local structure, peaks near -2, and reverses
            // outright on structureless high-dimensional data. See the StressMDS class docs.
            const X = swissroll(80);

            /** @param {number[][]} Y */
            const distance_correlation = (Y) => {
                /** @type {number[]} */ const input = [];
                /** @type {number[]} */ const output = [];
                for (let i = 0; i < X.length; ++i) {
                    for (let j = i + 1; j < X.length; ++j) {
                        input.push(Math.hypot(...X[i].map((v, c) => v - X[j][c])));
                        output.push(Math.hypot(...Y[i].map((v, c) => v - Y[j][c])));
                    }
                }
                const mean = (/** @type {number[]} */ z) => z.reduce((a, b) => a + b, 0) / z.length;
                const mi = mean(input);
                const mo = mean(output);
                let num = 0;
                let di = 0;
                let dout = 0;
                for (let i = 0; i < input.length; ++i) {
                    num += (input[i] - mi) * (output[i] - mo);
                    di += (input[i] - mi) ** 2;
                    dout += (output[i] - mo) ** 2;
                }
                return num / Math.sqrt(di * dout);
            };

            const [uniform, sammon, elastic] = [WEIGHTS_UNIFORM, WEIGHTS_SAMMON, WEIGHTS_ELASTIC].map((q) =>
                distance_correlation(new StressMDS(X, { weights: q, seed: 42 }).transform()),
            );
            expect(uniform).toBeGreaterThan(sammon);
            expect(sammon).toBeGreaterThan(elastic);
        });

        it("should expose the named exponents at their literature values", () => {
            expect(WEIGHTS_UNIFORM).toBe(0);
            expect(WEIGHTS_SAMMON).toBe(-1);
            expect(WEIGHTS_ELASTIC).toBe(-2);
        });

        it("should agree with SAMMON's objective at weights = -1", () => {
            // Sammon stress IS sigma_{-1} up to a constant factor, so StressMDS should reach a
            // comparable-or-lower value than SAMMON's own optimizer, from the same MDS start.
            const X = swissroll(80);
            const mine = weighted_stress(
                new StressMDS(X, { weights: WEIGHTS_SAMMON, seed: 42 }).transform(),
                X,
                -1,
            );
            const theirs = weighted_stress(new SAMMON(X, { seed: 42, init_DR: "MDS" }).transform(), X, -1);
            expect(mine).toBeLessThanOrEqual(theirs * 1.001);
        });

        it("should accept an explicit weight matrix", () => {
            const X = swissroll(30);
            const W = new Matrix(X.length, X.length, () => 1);
            const explicit = new StressMDS(X, { weights: W, seed: 42 }).transform();
            const uniform = new StressMDS(X, { weights: WEIGHTS_UNIFORM, seed: 42 }).transform();
            // w = 1 everywhere is exactly the q = 0 weighting.
            expect(explicit).toEqual(uniform);
        });

        it("should accept a plain nested array as the weight matrix", () => {
            // The Matrix path and the number[][] path go through different branches of _build_weights.
            const X = swissroll(30);
            const N = X.length;
            const nested = Array.from({ length: N }, (_, i) =>
                Array.from({ length: N }, (_, j) => (i === j ? 0 : 1)),
            );
            const asArray = new StressMDS(X, { weights: nested, seed: 42 }).transform();
            const asMatrix = new StressMDS(X, { weights: new Matrix(N, N, 1), seed: 42 }).transform();
            expect(asArray).toEqual(asMatrix);
        });

        it("should accept a weight function", () => {
            // `d ** -2` rather than `1 / (d * d)`: the two differ in the last bit, and the descent
            // amplifies that into a visible difference. Matching the exponent path's own expression
            // keeps this a test of routing rather than of floating point.
            const X = swissroll(30);
            const fn = new StressMDS(X, { weights: (d) => d ** -2, seed: 42 }).transform();
            const exponent = new StressMDS(X, { weights: WEIGHTS_ELASTIC, seed: 42 }).transform();
            expect(fn).toEqual(exponent);
        });

        it("should pass the pair indices to a weight function", () => {
            const X = swissroll(20);
            /** @type {Set<string>} */
            const seen = new Set();
            new StressMDS(X, {
                weights: (d, i, j) => {
                    seen.add(`${i},${j}`);
                    return d ** -1;
                },
                seed: 42,
            }).check_init();
            expect(seen.has("0,1")).toBe(true);
            expect(seen.has("19,18")).toBe(true);
            expect(seen.has("5,5")).toBe(false); // the diagonal is never asked about
        });

        it("should drop pairs whose weight is zero", () => {
            // The missing-data case: a zero weight must remove the pair from the objective, so those
            // pairs are free to land anywhere and the retained ones fit better than under w = 1.
            const X = swissroll(40);
            const N = X.length;
            const observed = (i, j) => (i + j) % 3 !== 0;
            const W = new Matrix(N, N, (i, j) => (i !== j && observed(i, j) ? 1 : 0));

            const partial = new StressMDS(X, { weights: W, seed: 42 }).transform();
            const full = new StressMDS(X, { weights: WEIGHTS_UNIFORM, seed: 42 }).transform();

            /** @param {number[][]} Y */
            const observed_stress = (Y) => {
                let total = 0;
                for (let i = 0; i < N; ++i)
                    for (let j = i + 1; j < N; ++j) {
                        if (!observed(i, j)) continue;
                        const d = Math.hypot(...X[i].map((v, c) => v - X[j][c]));
                        const r = Math.hypot(...Y[i].map((v, c) => v - Y[j][c])) - d;
                        total += r * r;
                    }
                return total;
            };

            expect(observed_stress(partial)).toBeLessThan(observed_stress(full));
        });

        it("should tolerate zero distances under a negative exponent", () => {
            // d = 0 gives 0^-2 = Infinity. Those pairs must be dropped, not given infinite stiffness.
            const X = swissroll(25);
            X.push([...X[0]], [...X[0]], [...X[3]]);
            const minfo = new StressMDS(X, { weights: WEIGHTS_ELASTIC, seed: 42 });
            expectValidValues(minfo.transform(), "StressMDS");
            expect(minfo.weights.entry(0, X.length - 3)).toBe(0);
        });

        it("should reject a mis-shaped weight matrix", () => {
            expect(() =>
                new StressMDS(generateTestData(20, 4), { weights: new Matrix(5, 5, 1) }).transform(),
            ).toThrow(/weights matrix/);
        });

        it("should reject an unusable weights value", () => {
            expect(() =>
                new StressMDS(generateTestData(20, 4), { weights: /** @type {any} */ ("nope") }).transform(),
            ).toThrow(/weights must be/);
        });
    });

    describe("the optimizer", () => {
        it("should decrease the stress monotonically", () => {
            const X = swissroll(40);
            const dr = new StressMDS(X, { seed: 42 });
            dr.check_init();
            const initial = dr.energy;

            let previous = Infinity;
            let steps = 0;
            for (const _ of dr.generator()) {
                expect(dr.energy).toBeLessThanOrEqual(previous);
                previous = dr.energy;
                steps++;
            }
            expect(steps).toBeGreaterThan(1);
            expect(dr.energy).toBeLessThan(initial);
        });

        it("should be insensitive to learning_rate", () => {
            // What the Jacobi preconditioning buys: the step is dimensionless, so a 20x change in
            // learning_rate must not meaningfully change where the descent lands.
            const X = swissroll(60);
            const values = [0.1, 0.5, 1, 2].map(
                (learning_rate) =>
                    new StressMDS(X, { learning_rate, seed: 42 }).transform() &&
                    new StressMDS(X, { learning_rate, seed: 42 }).energy,
            );
            const best = Math.min(...values);
            for (const value of values) expect(value / best).toBeLessThan(1.05);
        });

        it("should always yield at least the warm start", () => {
            const yielded = [...new StressMDS(generateTestData(20, 4), { seed: 42, iterations: 0 }).generator()];
            expect(yielded.length).toBeGreaterThanOrEqual(1);
            expect(yielded[0]).toHaveLength(20);
        });

        it("should support every init_DR", () => {
            // Only the PCA *rejection* was covered before, never a successful PCA start.
            const X = swissroll(40);
            for (const init_DR of /** @type {const} */ (["MDS", "PCA", "random"])) {
                const dr = new StressMDS(X, { init_DR, seed: 42 });
                const Y = dr.transform();
                expectValidValues(Y, `StressMDS(init_DR: ${init_DR})`);
                expect(Number.isFinite(dr.energy)).toBe(true);
            }
        });

        it("should reject PCA init on a precomputed matrix", () => {
            const D = new Matrix(10, 10, (i, j) => Math.abs(i - j));
            expect(() => new StressMDS(D, { metric: "precomputed", init_DR: "PCA" }).transform()).toThrow(
                /precomputed/,
            );
        });
    });

    describe("the WASM kernel", () => {
        it("should agree with the JS implementation", () => {
            const X = swissroll(80);

            setWasmEnabled(true);
            const accelerated = new StressMDS(X, { seed: 42 });
            const Y_wasm = accelerated.transform();

            setWasmEnabled(false);
            const plain = new StressMDS(X, { seed: 42 });
            const Y_js = plain.transform();

            for (let i = 0; i < Y_js.length; ++i) {
                for (let k = 0; k < Y_js[i].length; ++k) {
                    expect(Y_wasm[i][k]).toBeCloseTo(Y_js[i][k], 6);
                }
            }
            expect(accelerated.energy).toBeCloseTo(plain.energy, 4);
        });

        it("should agree with the JS implementation for a sparse weight matrix", () => {
            // The zero-weight skip is a branch in both paths; make sure they skip the same pairs.
            const X = swissroll(60);
            const N = X.length;
            const W = new Matrix(N, N, (i, j) => (i !== j && (i * j) % 4 !== 0 ? 1 / (i + j + 1) : 0));

            setWasmEnabled(true);
            const Y_wasm = new StressMDS(X, { weights: W, seed: 42 }).transform();
            setWasmEnabled(false);
            const Y_js = new StressMDS(X, { weights: W, seed: 42 }).transform();

            for (let i = 0; i < Y_js.length; ++i) {
                for (let k = 0; k < Y_js[i].length; ++k) {
                    expect(Y_wasm[i][k]).toBeCloseTo(Y_js[i][k], 6);
                }
            }
        });

        it("should release its buffers when a run ends", async () => {
            const { held_session_count } = await import("../../src/wasm/index.js");
            setWasmEnabled(true);
            new StressMDS(swissroll(40), { seed: 42 }).transform();
            expect(held_session_count()).toBe(0);
        });
    });

    describe("KKMDS preset", () => {
        it("should equal StressMDS at weights = -2", () => {
            const X = swissroll(40);
            expect(new KKMDS(X, { seed: 42 }).transform()).toEqual(
                new StressMDS(X, { weights: WEIGHTS_ELASTIC, seed: 42 }).transform(),
            );
        });

        it("should ignore an attempt to override weights", () => {
            const X = swissroll(40);
            const forced = new KKMDS(X, { .../** @type {any} */ ({ weights: 0 }), seed: 42 });
            forced.check_init();
            expect(forced.parameter("weights")).toBe(WEIGHTS_ELASTIC);
        });

        it("should work through its own statics", async () => {
            const X = generateTestData(20, 4);
            expect(KKMDS.transform(X, { seed: 42 })).toHaveLength(20);
            expect(await KKMDS.transform_async(X, { seed: 42 })).toHaveLength(20);
            const gen = KKMDS.generator(X, { seed: 42 });
            let last;
            for (const value of gen) last = value;
            expect(last).toHaveLength(20);
        });
    });
});
