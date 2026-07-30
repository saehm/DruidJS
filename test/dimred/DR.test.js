import { describe, expect, it } from "vitest";
import { DR } from "../../src/dimred/DR.js";
import { PCA } from "../../src/dimred/PCA.js";
import { TSNE } from "../../src/dimred/TSNE.js";
import { Matrix } from "../../src/matrix/Matrix.js";

/**
 * The `DR` base class carries contracts every method inherits — input-type round-tripping,
 * `parameter()`, the static convenience entry points, the generator protocol and WASM buffer
 * release. Each subclass test exercises its own algorithm, so the shared behaviour underneath was
 * only ever covered incidentally, and the parts a subclass overrides (`transform`, `generator`)
 * were never reached on the base at all.
 */

/** @param {number} n @param {number} d */
function grid(n, d) {
    return Array.from({ length: n }, (_, i) => Array.from({ length: d }, (_, j) => Math.sin(i * d + j) * 10 + i));
}

describe("DR (base class)", () => {
    describe("input type round-tripping", () => {
        const values = grid(20, 5);

        it("should return plain arrays for plain array input", () => {
            const Y = new PCA(values, { d: 2, seed: 1212 }).transform();
            expect(Array.isArray(Y)).toBe(true);
            expect(Array.isArray(Y[0])).toBe(true);
            expect(Y).toHaveLength(20);
            expect(Y[0]).toHaveLength(2);
        });

        it("should return Float64Arrays for Float64Array input", () => {
            const typed = values.map((row) => Float64Array.from(row));
            const Y = new PCA(typed, { d: 2, seed: 1212 }).transform();
            expect(Array.isArray(Y)).toBe(true);
            expect(Y[0]).toBeInstanceOf(Float64Array);
            expect(Y[0]).toHaveLength(2);
        });

        it("should return a Matrix for Matrix input", () => {
            const Y = new PCA(Matrix.from(values), { d: 2, seed: 1212 }).transform();
            expect(Y).toBeInstanceOf(Matrix);
            expect(/** @type {Matrix} */ (Y).shape).toEqual([20, 2]);
        });

        it("should produce the same numbers whichever container was used", () => {
            const from_array = new PCA(values, { d: 2, seed: 1212 }).transform();
            const from_matrix = /** @type {Matrix} */ (new PCA(Matrix.from(values), { d: 2, seed: 1212 }).transform());
            for (let i = 0; i < 20; ++i) {
                for (let j = 0; j < 2; ++j) {
                    expect(from_array[i][j]).toBeCloseTo(from_matrix.entry(i, j), 10);
                }
            }
        });

        it("should reject input that is neither an array nor a Matrix", () => {
            expect(() => new PCA(/** @type {any} */ ("not data"), { d: 2 })).toThrow("No valid type for X!");
            expect(() => new PCA(/** @type {any} */ (42), { d: 2 })).toThrow("No valid type for X!");
        });
    });

    describe("parameter()", () => {
        it("should return every parameter when called with no arguments", () => {
            const dr = new PCA(grid(10, 4), { d: 3, seed: 7 });
            const params = dr.parameter();
            expect(params.d).toBe(3);
            expect(params.seed).toBe(7);
        });

        it("should hand back a copy, so mutating it cannot reach the instance", () => {
            const dr = new PCA(grid(10, 4), { d: 3 });
            const params = dr.parameter();
            params.d = 99;
            expect(dr.parameter("d")).toBe(3);
        });

        it("should read a single parameter", () => {
            const dr = new PCA(grid(10, 4), { d: 3 });
            expect(dr.parameter("d")).toBe(3);
        });

        it("should set a parameter and return the instance for chaining", () => {
            const dr = new PCA(grid(10, 4), { d: 3 });
            expect(dr.parameter("d", 2)).toBe(dr);
            expect(dr.parameter("d")).toBe(2);
        });

        it("should re-run init after a parameter changes", () => {
            const data = grid(20, 5);
            const dr = new PCA(data, { d: 3, seed: 1212 });
            expect(/** @type {number[][]} */ (dr.transform())[0]).toHaveLength(3);

            // Setting a parameter clears `_is_initialized`, so the next transform must rebuild
            // rather than hand back the stale 3-column projection.
            dr.parameter("d", 2);
            expect(/** @type {number[][]} */ (dr.transform())[0]).toHaveLength(2);
        });

        it("should reject an unknown parameter name", () => {
            const dr = new PCA(grid(10, 4), { d: 2 });
            expect(() => dr.parameter(/** @type {any} */ ("nonexistent"))).toThrow("is not a valid parameter!");
            expect(() => dr.parameter(/** @type {any} */ ("nonexistent"), 1)).toThrow("is not a valid parameter!");
        });
    });

    describe("projection", () => {
        it("should throw when read before any transform", () => {
            const dr = new PCA(grid(10, 4), { d: 2 });
            expect(() => dr.projection).toThrow("The dataset is not transformed yet!");
        });

        it("should be readable again after transform", () => {
            const dr = new PCA(grid(10, 4), { d: 2, seed: 1212 });
            const Y = dr.transform();
            expect(dr.projection).toEqual(Y);
        });
    });

    describe("generator protocol", () => {
        it("should end with the same projection transform returns", () => {
            const data = grid(30, 5);
            const steps = [...new TSNE(data, { d: 2, perplexity: 5, seed: 1212 }).generator(10)];
            expect(steps.length).toBeGreaterThan(1);

            const direct = new TSNE(data, { d: 2, perplexity: 5, seed: 1212 }).transform(10);
            expect(steps.at(-1)).toEqual(direct);
        });

        it("should yield a projection of the right shape at every step", () => {
            for (const step of new TSNE(grid(30, 5), { d: 2, perplexity: 5, seed: 1212 }).generator(5)) {
                expect(step).toHaveLength(30);
                expect(step[0]).toHaveLength(2);
            }
        });
    });

    describe("async entry points", () => {
        it("should resolve to the same result as the synchronous transform", async () => {
            const data = grid(20, 5);
            const sync = new PCA(data, { d: 2, seed: 1212 }).transform();
            const async_result = await new PCA(data, { d: 2, seed: 1212 }).transform_async();
            expect(async_result).toEqual(sync);
        });

        it("should resolve the static transform_async", async () => {
            const data = grid(20, 5);
            const sync = PCA.transform(data, { d: 2, seed: 1212 });
            expect(await PCA.transform_async(data, { d: 2, seed: 1212 })).toEqual(sync);
        });
    });

    describe("base class fallbacks", () => {
        // `DR` is instantiable and its own `transform`/`generator` are the no-op identity the
        // subclasses replace. Nothing else reaches them.
        it("should expose the untransformed data through the base transform", () => {
            const dr = new DR(Matrix.from(grid(10, 4)), { seed: 1212 }, {});
            dr.Y = Matrix.from(grid(10, 2));
            expect(dr.transform()).toBeInstanceOf(Matrix);
        });

        it("should yield once from the base generator", () => {
            const dr = new DR(Matrix.from(grid(10, 4)), { seed: 1212 }, {});
            dr.Y = Matrix.from(grid(10, 2));
            const steps = [...dr.generator()];
            expect(steps).toHaveLength(1);
            expect(steps[0]).toBeInstanceOf(Matrix);
        });

        it("should return itself from check_init", () => {
            const dr = new DR(Matrix.from(grid(10, 4)), { seed: 1212 }, {});
            expect(dr.check_init()).toBe(dr);
        });
    });

    describe("WASM session release", () => {
        it("should hold no session keys on a method without an accelerated step", () => {
            expect(new PCA(grid(10, 4), { d: 2 })._wasm_session_keys).toEqual([]);
        });

        it("should name its sessions on a method with an accelerated step", () => {
            expect(new TSNE(grid(10, 4), { d: 2, perplexity: 3 })._wasm_session_keys.length).toBeGreaterThan(0);
        });

        it("should be safe to release when WASM never ran", () => {
            const dr = new PCA(grid(10, 4), { d: 2 });
            expect(() => dr._release_wasm()).not.toThrow();
            expect(() => dr._release_wasm()).not.toThrow();
        });

        it("should release through Symbol.dispose without disturbing the result", () => {
            const data = grid(30, 5);
            const dr = new TSNE(data, { d: 2, perplexity: 5, seed: 1212 });
            const Y = dr.transform(10);
            expect(() => dr[Symbol.dispose]()).not.toThrow();
            // Disposal frees WASM buffers, not the projection.
            expect(dr.projection).toEqual(Y);
        });

        it("should still project correctly after an abandoned generator is disposed", () => {
            const data = grid(30, 5);
            const partial = new TSNE(data, { d: 2, perplexity: 5, seed: 1212 });
            const it_ = partial.generator(10);
            it_.next();
            it_.next();
            partial[Symbol.dispose]();

            // A fresh run must not be affected by buffers released mid-flight elsewhere.
            const Y = new TSNE(data, { d: 2, perplexity: 5, seed: 1212 }).transform(10);
            expect(Y).toHaveLength(30);
            expect(Y.flat().every(Number.isFinite)).toBe(true);
        });
    });
});
