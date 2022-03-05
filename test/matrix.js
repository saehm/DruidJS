import * as druid from "./test_index.js";
//import "../dist/druid.js";
import * as assert from "assert";

describe("Matrix", () => {
    it("Create identity matrix", () => {
        let I = new druid.Matrix(3, 3, "I");
        assert.deepEqual(I.shape, [3, 3]);
        assert.equal(I.values.length, 3 * 3);
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) {
                assert.equal(I.entry(i, j), i == j ? 1 : 0);
            }
        }
    });

    it("Dot product", () => {
        let A = new druid.Matrix(5, 5, () => Math.random());
        let B = new druid.Matrix(5, 2, () => Math.random());

        assert.ok(A.dot(B));
        let C = A.dot(B);
        for (let i = 0; i < C.shape[0]; ++i) {
            for (let j = 0; j < C.shape[1]; ++j) {
                assert.ok(!Number.isNaN(C.entry(i, j)));
            }
        }
        assert.deepEqual(A.dot(B).shape, [5, 2]);
        assert.throws(() => B.dot(A), Error, "error thrown");
    });

    it("LU decomposition", () => {
        let A = new druid.Matrix(10, 10, () => Math.random());
        const { L, U } = druid.Matrix.LU(A);
        assert.deepEqual(L.shape, A.shape);
        assert.deepEqual(U.shape, A.shape);
        for (let row = 0; row < 10; ++row) {
            for (let col = 0; col < 10; ++col) {
                if (row < col) {
                    assert.equal(L.entry(row, col), 0);
                }
                if (row > col) {
                    assert.equal(U.entry(row, col), 0);
                }
            }
        }
    });
});

describe("norm", () => {
    it("norm", () => {
        const v = druid.normalize(Float64Array.from({length: 100}, () => Math.random() * 100 - 50));
        assert.equal(Math.abs(druid.norm(v) - 1) < 1e-12, 1)
        
    })
})