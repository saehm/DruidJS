//import "../dist/druid.js";
import * as druid from "./test_index.js";
import * as assert from "assert";

const eps = 0.00000001;
describe("eig", () => {
    const N = 20;
    const R = new druid.Randomizer(1212);
    let M = new druid.Matrix(N, N, () => (R.random - .5)  * 5);
    M = M.dot(M.T);

    it("qr", () => {
        assert.ok(druid.qr(M));
        assert.ok(druid.qr_householder(M));
        const { Q: QM, R: RM } = druid.qr(M);
        checkDecomposition(M, QM, RM);
        const { Q: QMh, R: RMh } = druid.qr_householder(M);
        checkDecomposition(M, QMh, RMh);

        const A = druid.Matrix.from([[12, -51, 4], [6, 167, -68], [-4, 24, -41]]);
        const { Q, R } = druid.qr(A);
        approxEqual(Q.values, Float64Array.from([
          6/7, -69/175, -58/175,
          3/7, 158/175, 6/175,
          -2/7, 6/35, -33/35
        ]));
        approxEqual(R.values, Float64Array.from([14, 21, -14, 0, 175, -70, 0, 0, 35]));

        const B = druid.Matrix.from([
          [1, -1,  4],
          [1,  4, -2],
          [1,  4,  2],
          [1,  -1, 0]
        ]);
        const { Q: QB, R: RB } = druid.qr(B);
        approxEqual(QB.values, Float64Array.from([
            0.5, -0.5,  0.5,
            0.5,  0.5, -0.5,
            0.5,  0.5,  0.5,
            0.5, -0.5, -0.5
        ]));
        approxEqual(RB.values, Float64Array.from([2, 3, 2, 0, 5, -2, 0, 0, 4]));
    });
    it("simultanious poweriteration", () => {
        assert.ok(druid.simultaneous_poweriteration(M, 2, {max_iterations: 100, seed: 1212, qr: druid.qr}));
        assert.ok(druid.simultaneous_poweriteration(M, 2, {max_iterations: 100, seed: 1212, qr: druid.qr_householder}));
        let eigs = druid.simultaneous_poweriteration(M, 2, {max_iterations: 100, seed: 1212, qr: druid.qr});
        assert.equal(eigs.eigenvectors.length, 2)
        assert.equal(eigs.eigenvectors[0].length, N)

        eigs = druid.simultaneous_poweriteration(M, 6, {max_iterations: 100, seed: 1212, qr: druid.qr});
        assert.equal(eigs.eigenvectors.length, 6)
        assert.equal(eigs.eigenvectors[0].length, N)
    }).timeout(10000);
});

function checkDecomposition(A, Q, R) {
    approxEqual(Q.dot(R).values, A.values);

    // R is upper triangular
    const [rows, cols] = A.shape;
    for (let i = 0; i < rows; ++i) {
        for (let j = 0; j < i && j < cols; ++j) {
            assert.ok(Math.abs(R.entry(i, j)) < eps);
        }
    }

    // All elements on leading diagonal of R are positive
    // for (let i = 0; i < Math.min(rows, cols); i++) {
    //     assert.ok(R.entry(i, i) >= 0);
    // }
}

function approxEqual(a, b) {
  const N = a.length;
  assert.ok(a.length === b.length);
  for (let i = 0; i < N; ++i) {
    assert.ok(Math.abs(a[i] - b[i]) < eps);
  }
}