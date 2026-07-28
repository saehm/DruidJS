import { describe, expect, it } from "vitest";
import { Matrix } from "../../src/matrix/index.js";
import { setWasmEnabled } from "../../src/wasm/index.js";

/**
 * `A · x` computed directly from the definition, as the reference.
 *
 * @param {number[][]} A
 * @param {number[]} x
 * @returns {number[]}
 */
function mat_vec(A, x) {
    return A.map((row) => row.reduce((sum, a, j) => sum + a * x[j], 0));
}

describe("Matrix vector products", () => {
    const A2 = [
        [1, 2],
        [3, 4],
    ];

    it("dot(array) computes A · x", () => {
        const C = Matrix.from(A2).dot([10, 20]);
        expect(C.shape).toEqual([2, 1]);
        expect(C.to2dArray().map((r) => [...r])).toEqual([[50], [110]]);
    });

    it("dot(array) agrees with dot(Matrix)", () => {
        const A = Matrix.from(A2);
        const as_column = Matrix.from([[10], [20]]);
        expect(
            A.dot([10, 20])
                .to2dArray()
                .map((r) => [...r]),
        ).toEqual(
            A.dot(as_column)
                .to2dArray()
                .map((r) => [...r]),
        );
    });

    it("transDot(array) computes Aᵀ · x", () => {
        // Aᵀ = [[1, 3], [2, 4]] so Aᵀ · [10, 20] = [1*10+3*20, 2*10+4*20]
        const C = Matrix.from(A2).transDot([10, 20]);
        expect(C.shape).toEqual([2, 1]);
        expect(C.to2dArray().map((r) => [...r])).toEqual([[70], [100]]);
    });

    it("dotTrans(array) computes A · xᵀ", () => {
        expect(
            Matrix.from(A2)
                .dotTrans([10, 20])
                .to2dArray()
                .map((r) => [...r]),
        ).toEqual([[50], [110]]);
    });

    it("rejects a vector whose length does not match", () => {
        const A = Matrix.from([[1, 2, 3]]); // 1 ⨯ 3
        expect(() => A.dot([1, 2])).toThrow();
        expect(() => A.transDot([1, 2])).toThrow(); // needs `rows` = 1 entries
    });

    it("handles a non-square matrix", () => {
        const A = [
            [1, 2, 3],
            [4, 5, 6],
        ];
        expect(
            Matrix.from(A)
                .dot([1, 1, 1])
                .to2dArray()
                .map((r) => [...r]),
        ).toEqual([[6], [15]]);
        expect(
            Matrix.from(A)
                .transDot([1, 1])
                .to2dArray()
                .map((r) => [...r]),
        ).toEqual([[5], [7], [9]]);
    });

    it("accepts a Float64Array", () => {
        const C = Matrix.from(A2).dot(new Float64Array([10, 20]));
        expect(C.to2dArray().map((r) => [...r])).toEqual([[50], [110]]);
    });

    // Large enough to cross WASM_MIN_MATMUL_OPS, so both paths are exercised.
    it("WASM and JS paths agree above the kernel threshold", () => {
        const rows = 40;
        const cols = 40;
        const A = Array.from({ length: rows }, (_, i) =>
            Array.from({ length: cols }, (_, j) => Math.sin(i * cols + j)),
        );
        const x = Array.from({ length: cols }, (_, j) => Math.cos(j));
        const expected = mat_vec(A, x);

        setWasmEnabled(true);
        const with_wasm = Matrix.from(A)
            .dot(x)
            .to2dArray()
            .map((r) => r[0]);
        setWasmEnabled(false);
        const without_wasm = Matrix.from(A)
            .dot(x)
            .to2dArray()
            .map((r) => r[0]);
        setWasmEnabled(true);

        for (let i = 0; i < rows; ++i) {
            expect(with_wasm[i]).toBeCloseTo(expected[i], 10);
            expect(without_wasm[i]).toBeCloseTo(expected[i], 10);
        }
    });
});
