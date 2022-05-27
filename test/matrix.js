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

        const D = new druid.Matrix(2, 3);
        D._data = Float64Array.from([0, 1, 2, 3, 4, 5]);

        assert.ok(D.dot(D.T));
        assert.ok(D.dotTrans(D));
        const D_dot_DT = Float64Array.from([5, 14, 14, 50]);
        assert.deepEqual(D.dot(D.T).values, D_dot_DT);
        assert.deepEqual(D.dotTrans(D).values, D_dot_DT);

        assert.ok(D.T.dot(D));
        assert.ok(D.transDot(D));
        const DT_dot_D = Float64Array.from([9, 12, 15, 12, 17, 22, 15, 22, 29]);
        assert.deepEqual(D.T.dot(D).values, DT_dot_D);
        assert.deepEqual(D.transDot(D).values, DT_dot_D);
    });

    it("Matrix inversion", () => {
      const A = new druid.Matrix(2, 2);
      A._data = Float64Array.from([2, 3, 4, 7]);
      assert.ok(A.inverse());
      assert.deepEqual(A.inverse().values, Float64Array.from([3.5, -1.5, -2, 1]));
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