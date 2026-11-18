import { describe, expect, it } from "vitest";
import { SMACOF } from "../../src/dimred/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { generateTestData } from "../utils/data-generators.js";
import { expectValidValues } from "../utils/helpers.js";

describe("SMACOF", () => {
    it("should reduce dimensionality to specified d", () => {
        const data = generateTestData(15, 4);
        const smacof = new SMACOF(data, { d: 2, seed: 42 });
        const result = smacof.transform();

        expect(result).toHaveLength(15);
        expect(result[0]).toHaveLength(2);
    });

    it("should produce valid values", { timeout: 10000 }, () => {
        const data = generateTestData(15, 5);
        const smacof = new SMACOF(data, { d: 2, seed: 42 });
        const result = smacof.transform();
        expectValidValues(result, "SMACOF");
    });

    it("should work with static transform", () => {
        const data = generateTestData(10, 5);
        const result = SMACOF.transform(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with static generator", () => {
        const data = generateTestData(10, 5);
        const gen = SMACOF.generator(data, { d: 2 });
        const result = gen.next().value;
        expect(result).toHaveLength(10);
    });

    it("should work with static transform_async", async () => {
        const data = generateTestData(10, 5);
        const result = await SMACOF.transform_async(data, { d: 2 });
        expect(result).toHaveLength(10);
    });

    it("should work with precomputed metric", () => {
        const data = generateTestData(10, 5);
        const dists = Matrix.from(generateTestData(10, 10)); // dummy dist matrix
        const smacof = new SMACOF(dists, { metric: "precomputed", d: 2 });
        const result = smacof.transform();
        expect(result.shape[0]).toBe(10);
    });

    describe("start and continuation", () => {
        const data = generateTestData(20, 4);
        /** @param {Matrix} A @param {Matrix} B */
        const max_diff = (A, B) => Math.max(...A.values.map((v, i) => Math.abs(v - B.values[i])));

        it("honors the iteration count passed to generator()", () => {
            const smacof = new SMACOF(data, { seed: 42, epsilon: 0 });
            expect([...smacof.generator(5)]).toHaveLength(5);
            expect([...new SMACOF(data, { seed: 42, epsilon: 0, iterations: 3 }).generator()]).toHaveLength(3);
        });

        it("starts from init_Y", () => {
            const init_Y = new Matrix(20, 2, (i, j) => (j === 0 ? i : (i * 7) % 5));
            const smacof = new SMACOF(Matrix.from(data), { seed: 42, init_Y });
            const first = /** @type {Matrix} */ (smacof.generator(1).next().value);
            // one Guttman step from init_Y, computed by hand
            const D = new Matrix(20, 20, (i, j) => Math.hypot(...data[i].map((v, k) => v - data[j][k])));
            const step = new Matrix(20, 2, (i, c) => {
                let s = 0;
                for (let j = 0; j < 20; ++j) {
                    if (i === j) continue;
                    const dz = Math.hypot(
                        init_Y.entry(i, 0) - init_Y.entry(j, 0),
                        init_Y.entry(i, 1) - init_Y.entry(j, 1),
                    );
                    const b = dz > 1e-12 ? D.entry(i, j) / dz : 0;
                    s += b * (init_Y.entry(i, c) - init_Y.entry(j, c));
                }
                return s / 20;
            });
            expect(max_diff(first, step)).toBeLessThan(1e-10);
            expect(init_Y.entry(3, 0)).toBe(3); // not modified
        });

        it("throws for an init_Y of the wrong shape", () => {
            const smacof = new SMACOF(data, { init_Y: new Matrix(20, 3, 0) });
            expect(() => smacof.transform()).toThrow(/init_Y needs shape 20 × 2/);
        });

        it("continues from the current layout when called again", () => {
            const X = Matrix.from(data);
            const once = new SMACOF(X, { seed: 42, epsilon: 0 }).transform(40);
            const twice = new SMACOF(X, { seed: 42, epsilon: 0 });
            twice.transform(15);
            const result = twice.transform(25);
            expect(max_diff(result, once)).toBeLessThan(1e-12);
        });
    });
});
