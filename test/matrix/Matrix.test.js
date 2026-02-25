import { describe, expect, test } from "vitest";
import { Matrix } from "../../src/matrix/index.js";

describe("Matrix", () => {
    describe("Constructor", () => {
        test("creates empty matrix with no value", () => {
            const m = new Matrix(2, 3);
            expect(m.shape).toEqual([2, 3]);
            expect(m.values.length).toBe(6);
            expect(m.values.every((v) => v === 0)).toBe(true);
        });

        test("creates matrix with number value", () => {
            const m = new Matrix(2, 2, 5);
            expect(m.values).toEqual(new Float64Array([5, 5, 5, 5]));
        });

        test("creates matrix with function", () => {
            const m = new Matrix(2, 3, (i, j) => i * 10 + j);
            expect(m.entry(0, 0)).toBe(0);
            expect(m.entry(0, 1)).toBe(1);
            expect(m.entry(1, 0)).toBe(10);
            expect(m.entry(1, 2)).toBe(12);
        });

        test("creates zeros matrix", () => {
            const m = new Matrix(2, 2, "zeros");
            expect(m.values.every((v) => v === 0)).toBe(true);
        });

        test("creates identity matrix", () => {
            const m = new Matrix(3, 3, "I");
            expect(m.entry(0, 0)).toBe(1);
            expect(m.entry(0, 1)).toBe(0);
            expect(m.entry(1, 1)).toBe(1);
            expect(m.entry(2, 2)).toBe(1);
        });

        test("creates center matrix", () => {
            const m = new Matrix(3, 3, "center");
            expect(m.entry(0, 0)).toBeCloseTo(2 / 3, 10);
            expect(m.entry(0, 1)).toBeCloseTo(-1 / 3, 10);
            expect(m.entry(1, 0)).toBeCloseTo(-1 / 3, 10);
        });

        test("creates matrix from 2D array", () => {
            const m = new Matrix(2, 2, [
                [1, 2],
                [3, 4],
            ]);
            expect(m.entry(0, 0)).toBe(1);
            expect(m.entry(0, 1)).toBe(2);
            expect(m.entry(1, 0)).toBe(3);
            expect(m.entry(1, 1)).toBe(4);
        });
    });

    describe("Static factory methods", () => {
        test("Matrix.from creates from 2D array", () => {
            const m = Matrix.from([
                [1, 2, 3],
                [4, 5, 6],
            ]);
            expect(m.shape).toEqual([2, 3]);
            expect(m.entry(0, 2)).toBe(3);
            expect(m.entry(1, 0)).toBe(4);
        });

        test("Matrix.from clones Matrix", () => {
            const original = new Matrix(2, 2, [
                [1, 2],
                [3, 4],
            ]);
            const clone = Matrix.from(original);
            expect(clone.shape).toEqual(original.shape);
            expect(clone._data).not.toBe(original._data);
            expect(clone.values).toEqual(original.values);
        });

        test("Matrix.from throws on invalid input", () => {
            expect(() => Matrix.from("invalid")).toThrow();
        });

        test("Matrix.from_diag creates diagonal matrix", () => {
            const m = Matrix.from_diag([1, 2, 3]);
            expect(m.shape).toEqual([3, 3]);
            expect(m.entry(0, 0)).toBe(1);
            expect(m.entry(1, 1)).toBe(2);
            expect(m.entry(2, 2)).toBe(3);
            expect(m.entry(0, 1)).toBe(0);
        });

        test("Matrix.from_vector creates column vector", () => {
            const m = Matrix.from_vector([1, 2, 3], "col");
            expect(m.shape).toEqual([3, 1]);
            expect(m.entry(0, 0)).toBe(1);
            expect(m.entry(2, 0)).toBe(3);
        });

        test("Matrix.from_vector creates row vector", () => {
            const m = Matrix.from_vector([1, 2, 3], "row");
            expect(m.shape).toEqual([1, 3]);
            expect(m.entry(0, 0)).toBe(1);
            expect(m.entry(0, 2)).toBe(3);
        });
    });

    describe("Row and column access", () => {
        test("row returns Float64Array", () => {
            const m = Matrix.from([
                [1, 2, 3],
                [4, 5, 6],
            ]);
            const row = m.row(1);
            expect(row).toBeInstanceOf(Float64Array);
            expect(row).toEqual(new Float64Array([4, 5, 6]));
        });

        test("iterate_rows yields all rows", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            const rows = [...m.iterate_rows()];
            expect(rows.length).toBe(2);
            expect(rows[0]).toEqual(new Float64Array([1, 2]));
            expect(rows[1]).toEqual(new Float64Array([3, 4]));
        });

        test("Symbol.iterator works", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            const rows = [];
            for (const row of m) {
                rows.push(row);
            }
            expect(rows.length).toBe(2);
        });

        test("col returns Float64Array", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
                [5, 6],
            ]);
            const col = m.col(1);
            expect(col).toBeInstanceOf(Float64Array);
            expect(col).toEqual(new Float64Array([2, 4, 6]));
        });

        test("set_row with array", () => {
            const m = new Matrix(2, 2);
            m.set_row(0, [1, 2]);
            expect(m.row(0)).toEqual(new Float64Array([1, 2]));
        });

        test("swap_rows exchanges rows", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            m.swap_rows(0, 1);
            expect(m.row(0)).toEqual(new Float64Array([3, 4]));
            expect(m.row(1)).toEqual(new Float64Array([1, 2]));
        });
    });

    describe("Entry operations", () => {
        test("entry gets value", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            expect(m.entry(0, 0)).toBe(1);
        });

        test("set_entry sets value", () => {
            const m = new Matrix(2, 2);
            m.set_entry(0, 1, 5);
            expect(m.entry(0, 1)).toBe(5);
        });

        test("add_entry adds to value", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            m.add_entry(0, 0, 10);
            expect(m.entry(0, 0)).toBe(11);
        });

        test("sub_entry subtracts from value", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            m.sub_entry(0, 0, 1);
            expect(m.entry(0, 0)).toBe(0);
        });
    });

    describe("Matrix operations", () => {
        test("transpose swaps dimensions", () => {
            const m = Matrix.from([
                [1, 2, 3],
                [4, 5, 6],
            ]);
            const t = m.transpose();
            expect(t.shape).toEqual([3, 2]);
            expect(t.entry(0, 1)).toBe(4);
            expect(t.entry(2, 0)).toBe(3);
        });

        test("dot multiplies matrices", () => {
            const a = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            const b = Matrix.from([
                [5, 6],
                [7, 8],
            ]);
            const c = a.dot(b);
            expect(c.shape).toEqual([2, 2]);
            expect(c.entry(0, 0)).toBe(19);
            expect(c.entry(1, 1)).toBe(50);
        });

        test("outer computes outer product", () => {
            const a = Matrix.from_vector([1, 2, 3], "col");
            const b = Matrix.from_vector([4, 5, 6], "col");
            const c = a.outer(b);
            expect(c.shape).toEqual([3, 3]);
            expect(c.entry(0, 0)).toBe(4);
            expect(c.entry(2, 2)).toBe(18);
        });
    });

    describe("Concatenation", () => {
        test("concat horizontal", () => {
            const a = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            const b = Matrix.from([
                [5, 6],
                [7, 8],
            ]);
            const c = a.concat(b, "horizontal");
            expect(c.shape).toEqual([2, 4]);
        });

        test("concat vertical", () => {
            const a = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            const b = Matrix.from([
                [5, 6],
                [7, 8],
            ]);
            const c = a.concat(b, "vertical");
            expect(c.shape).toEqual([4, 2]);
        });
    });

    describe("Block operations", () => {
        test("get_block extracts submatrix", () => {
            const m = Matrix.from([
                [1, 2, 3],
                [4, 5, 6],
                [7, 8, 9],
            ]);
            const block = m.get_block(0, 0, 2, 2);
            expect(block.shape).toEqual([2, 2]);
            expect(block.entry(0, 0)).toBe(1);
            expect(block.entry(1, 1)).toBe(5);
        });

        test("set_block inserts matrix", () => {
            const m = new Matrix(3, 3);
            const block = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            m.set_block(1, 1, block);
            expect(m.entry(1, 1)).toBe(1);
            expect(m.entry(2, 2)).toBe(4);
        });
    });

    describe("Arithmetic operations", () => {
        test("mult with scalar", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            const result = m.mult(2);
            expect(result.values).toEqual(new Float64Array([2, 4, 6, 8]));
        });

        test("add with scalar", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            const result = m.add(10);
            expect(result.values).toEqual(new Float64Array([11, 12, 13, 14]));
        });

        test("sub with scalar", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            const result = m.sub(1);
            expect(result.values).toEqual(new Float64Array([0, 1, 2, 3]));
        });
    });

    describe("Properties and utilities", () => {
        test("clone creates independent copy", () => {
            const a = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            const b = a.clone();
            expect(b._data).not.toBe(a._data);
            expect(b.values).toEqual(a.values);
        });

        test("sum returns total", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            expect(m.sum()).toBe(10);
        });

        test("mean returns average", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            expect(m.mean()).toBe(2.5);
        });
    });

    describe("Inverse", () => {
        test("inverse of 2x2 matrix", () => {
            const m = Matrix.from([
                [4, 7],
                [2, 6],
            ]);
            const inv = m.inverse();
            const product = m.dot(inv);
            expect(product.entry(0, 0)).toBeCloseTo(1, 10);
            expect(product.entry(1, 1)).toBeCloseTo(1, 10);
        });
    });

    describe("Static methods", () => {
        test("Matrix.det computes determinant", () => {
            const m = Matrix.from([
                [4, 7],
                [2, 6],
            ]);
            expect(Matrix.det(m)).toBeCloseTo(10, 10);
        });

        test("Matrix.solve solves Ax = b", () => {
            const A = Matrix.from([
                [2, 5],
                [1, 3],
            ]);
            const b = Matrix.from([[12], [7]]);
            const x = Matrix.solve(A, b);
            const Ax = A.dot(x);
            expect(Ax.entry(0, 0)).toBeCloseTo(12, 6);
            expect(Ax.entry(1, 0)).toBeCloseTo(7, 6);
        });
    });

    describe("Additional Matrix Methods", () => {
        test("diag returns diagonal elements", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            expect(m.diag()).toEqual(new Float64Array([1, 4]));
        });

        test("meanRows and meanCols", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
            ]);
            expect(m.meanRows()).toEqual(new Float64Array([1.5, 3.5]));
            expect(m.meanCols()).toEqual(new Float64Array([2, 3]));
        });

        test("LU decomposition", () => {
            const m = Matrix.from([
                [4, 3],
                [6, 3],
            ]);
            const { L, U } = Matrix.LU(m);
            const LU = L.dot(U);
            expect(LU.entry(0, 0)).toBeCloseTo(4);
            expect(LU.entry(1, 1)).toBeCloseTo(3);
        });

        test("SVD decomposition", () => {
            const m = Matrix.from([
                [1, 2],
                [3, 4],
                [5, 6],
            ]);
            const { U, Sigma, V } = Matrix.SVD(m, 2);
            expect(Sigma.length).toBe(2);
            expect(V.length).toBe(2);
        });

        test("solve_CG solves Ax = b", () => {
            const A = Matrix.from([
                [4, 1],
                [1, 3],
            ]);
            const b = Matrix.from([[1], [2]]);
            const x = Matrix.solve_CG(A, b);
            const Ax = A.dot(x);
            expect(Ax.entry(0, 0)).toBeCloseTo(1, 2);
            expect(Ax.entry(1, 0)).toBeCloseTo(2, 2);
        });

        test("gather extracts specific rows and cols", () => {
            const m = Matrix.from([
                [1, 2, 3],
                [4, 5, 6],
                [7, 8, 9],
            ]);
            const result = m.gather([0, 2], [0, 2]);
            expect(result.shape).toEqual([2, 2]);
            expect(result.entry(0, 0)).toBe(1);
            expect(result.entry(1, 1)).toBe(9);
        });

        test("concat diag", () => {
            const a = Matrix.from([
                [1, 1],
                [1, 1],
            ]);
            const b = Matrix.from([
                [2, 2],
                [2, 2],
            ]);
            const c = a.concat(b, "diag");
            expect(c.shape).toEqual([4, 4]);
            expect(c.entry(0, 0)).toBe(1);
            expect(c.entry(2, 2)).toBe(2);
            expect(c.entry(0, 2)).toBe(0);
        });
    });

    describe("Static math methods", () => {
        test("LU throws on 0 diagonal", () => {
            const m = Matrix.from([
                [0, 1],
                [1, 0],
            ]);
            expect(() => Matrix.LU(m)).toThrow("L's diagonal not supposed to be 0!");
        });

        test("det 3x3", () => {
            const m = Matrix.from([
                [1, 2, 3],
                [4, 5, 6],
                [7, 8, 9],
            ]);
            expect(Matrix.det(m)).toBeCloseTo(0, 10);
        });

        test("det 4x4", () => {
            // An upper triangular 4x4 matrix avoids 0 pivot
            const m = Matrix.from([
                [1, 2, 3, 4],
                [0, 5, 6, 7],
                [0, 0, 8, 9],
                [0, 0, 0, 10],
            ]);
            expect(Matrix.det(m)).toBeCloseTo(400, 10);
        });
    });

    describe("Static validation methods", () => {
        test("Matrix.is2dArray validation", () => {
            expect(
                Matrix.is2dArray([
                    [1, 2],
                    [3, 4],
                ]),
            ).toBe(true);
            expect(Matrix.is2dArray([new Float64Array([1, 2]), new Float64Array([3, 4])])).toBe(true);
            expect(Matrix.is2dArray([])).toBe(false);
            expect(Matrix.is2dArray([[1, 2], [3]])).toBe(false);
            expect(Matrix.is2dArray([[1, 2], "invalid"])).toBe(false);
        });
    });
});
