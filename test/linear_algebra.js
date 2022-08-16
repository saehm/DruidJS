//import "../dist/druid.js";
import * as druid from "./test_index.js";
import * as assert from "assert";

const eps = 0.000001;
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
        // const { Q: QMh, R: RMh } = druid.qr_householder(M);
        // checkDecomposition(M, QMh, RMh);

        const M1 = druid.Matrix.from([[15, 42], [20, 81]]);
        const { Q: QM1, R: RM1 } = druid.qr(M1);
        approxEqual(QM1.values, Float64Array.from([0.6, -0.8, 0.8, 0.6]));
        approxEqual(RM1.values, Float64Array.from([25, 90, 0, 15]));
        checkDecomposition(M1, QM1, RM1);

        const M2 = druid.Matrix.from([[12, -51, 4], [6, 167, -68], [-4, 24, -41]]);
        const { Q: QM2, R: RM2 } = druid.qr(M2);
        approxEqual(QM2.values, Float64Array.from([
          6/7, -69/175, -58/175,
          3/7, 158/175, 6/175,
          -2/7, 6/35, -33/35
        ]));
        approxEqual(RM2.values, Float64Array.from([14, 21, -14, 0, 175, -70, 0, 0, 35]));
        checkDecomposition(M2, QM2, RM2);

        const M3 = druid.Matrix.from([
          [1, -1,  4],
          [1,  4, -2],
          [1,  4,  2],
          [1,  -1, 0]
        ]);
        const { Q: QM3, R: RM3 } = druid.qr(M3);
        approxEqual(QM3.values, Float64Array.from([
            0.5, -0.5,  0.5,
            0.5,  0.5, -0.5,
            0.5,  0.5,  0.5,
            0.5, -0.5, -0.5
        ]));
        approxEqual(RM3.values, Float64Array.from([2, 3, 2, 0, 5, -2, 0, 0, 4]));
        checkDecomposition(M3, QM3, RM3);

        const M4 = druid.Matrix.from([
          [7.507, 9.868, 5.057],
          [4.482, 2.536, 9.744],
          [6.527, 1.094, 3.321]
        ]);
        const { Q: QM4, R: RM4 } = druid.qr(M4);
        checkDecomposition(M4, QM4, RM4);
    });
    it("simultanious poweriteration", () => {
        assert.ok(druid.simultaneous_poweriteration(M, 2, {max_iterations: 100, seed: 1212, qr: druid.qr}));
        assert.ok(druid.simultaneous_poweriteration(M, 2, {max_iterations: 100, seed: 1212, qr: druid.qr_householder}));
        let eigs = druid.simultaneous_poweriteration(M, 2, {max_iterations: 100, seed: 1212, qr: druid.qr});
        assert.equal(eigs.eigenvectors.length, 2)
        assert.equal(eigs.eigenvectors[0].length, N)
        checkEigs(M, eigs.eigenvalues, eigs.eigenvectors);

        eigs = druid.simultaneous_poweriteration(M, 6, {max_iterations: 100, seed: 1212, qr: druid.qr});
        assert.equal(eigs.eigenvectors.length, 6)
        assert.equal(eigs.eigenvectors[0].length, N)
        checkEigs(M, eigs.eigenvalues, eigs.eigenvectors);

        const A = druid.Matrix.from([[1, 0.1], [0.1, 1]]);
        const { eigenvalues: A_val } = druid.simultaneous_poweriteration(A, 2);
        approxEqual(A_val, Float64Array.from([1.1, 0.9]));

        // const B = druid.Matrix.from([[3, -1, -1], [-12, 0, 5], [4, -2, -1]]);
        const B = druid.Matrix.from([[13, -4, 2], [-4, 11, -2], [2, -2, 8]]);
        const { eigenvalues: B_val, eigenvectors: B_vec } = druid.simultaneous_poweriteration(B, 3);
        approxEqual(B_val, Float64Array.from([17, 8, 7]));
        checkEigs(B, B_val, B_vec);
    }).timeout(10000);
});

function checkDecomposition(A, Q, R) {
    approxEqual(Q.dot(R).values, A.values);

    // R is upper triangular
    const [rows, cols] = A.shape;
    for (let i = 0; i < cols; ++i) {
        for (let j = 0; j < i; ++j) {
            assert.ok(Math.abs(R.entry(i, j)) < eps);
        }
    }

    // All elements on leading diagonal of R are positive
    for (let i = 0; i < Math.min(rows, cols); ++i) {
        assert.ok(R.entry(i, i) >= 0);
    }
}

function checkEigs(A, vals, vecs) {
  vals.forEach((val, i) => {
    const proj1 = A.dot(druid.Matrix.from(vecs[i], 'col'));
    const proj2 = vecs[i].map(d => d * val);
    approxEqual(proj1.values, proj2, 0.005);
  });
}

function approxEqual(a, b, e = eps) {
  // assert.deepEqual(a, b);
  const N = a.length;
  assert.ok(a.length === b.length);
  for (let i = 0; i < N; ++i) {
    assert.ok(Math.abs(a[i] - b[i]) < e, a[i] + ' != ' + b[i]);
  }
}