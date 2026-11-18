import { describe, expect, it } from "vitest";
import { inner_product, qr, qr_householder, simultaneous_poweriteration } from "../../src/linear_algebra/index.js";
import { Matrix } from "../../src/matrix/index.js";
import { Randomizer } from "../../src/util/index.js";

describe("inner_product", () => {
    it("should compute inner product", () => {
        expect(inner_product([1, 2, 3], [4, 5, 6])).toBe(32);
    });

    it("should return 0 for orthogonal vectors", () => {
        expect(inner_product([1, 0], [0, 1])).toBe(0);
    });
});

describe("qr", () => {
    it("should decompose 2x2 matrix", () => {
        const A = Matrix.from([
            [1, 1],
            [0, 1],
        ]);
        const { Q, R } = qr(A);
        const QR = Q.dot(R);
        expect(QR.entry(0, 0)).toBeCloseTo(1, 5);
        expect(R.entry(1, 0)).toBeCloseTo(0, 5);
    });
});

describe("qr_householder", () => {
    it("should decompose 2x2 matrix", () => {
        const A = Matrix.from([
            [1, 1],
            [0, 1],
        ]);
        const { Q, R } = qr_householder(A);
        const QR = Q.dot(R);
        expect(QR.entry(0, 0)).toBeCloseTo(1, 5);
        expect(R.entry(1, 0)).toBeCloseTo(0, 5);
    });
});

describe("simultaneous_poweriteration", () => {
    it("should compute eigenvalues and eigenvectors", () => {
        const A = Matrix.from([
            [2, 1],
            [1, 2],
        ]);
        const { eigenvalues, eigenvectors } = simultaneous_poweriteration(A, 2, { seed: 42 });
        expect(eigenvalues[0]).toBeCloseTo(3, 3);
        expect(eigenvalues[1]).toBeCloseTo(1, 3);
    });

    it("should keep the sign of negative eigenvalues", () => {
        const A = Matrix.from([
            [5, 0, 0],
            [0, -3, 0],
            [0, 0, 1],
        ]);
        const { eigenvalues } = simultaneous_poweriteration(A, 2, { seed: 42 });
        expect(eigenvalues[0]).toBeCloseTo(5, 6);
        expect(eigenvalues[1]).toBeCloseTo(-3, 6);
    });

    it("should stop early when an eigenvalue is negative", () => {
        // The eigenvector of a negative eigenvalue flips its sign in every step; convergence must
        // still be detected instead of running into max_iterations.
        /**
         * @param {number} d
         * @returns {number}
         */
        const steps = (d) => {
            let calls = 0;
            /** @type {typeof qr} */
            const counting_qr = (M) => {
                ++calls;
                return qr(M);
            };
            const A = Matrix.from([
                [5, 0, 0],
                [0, d, 0],
                [0, 0, 1],
            ]);
            simultaneous_poweriteration(A, 2, { seed: 42, qr: counting_qr, max_iterations: 1000 });
            return calls - 1; // the first call orthogonalizes the random start
        };
        expect(steps(-3)).toBeLessThan(100);
        expect(steps(-3)).toBeLessThan(2 * steps(3));
    });
});

describe("simultaneous_poweriteration with smallest: true", () => {
    /**
     * ‖A v − λ v‖ for an eigenpair.
     *
     * @param {Matrix} A
     * @param {Float64Array} v
     * @param {number} lambda
     */
    const residual = (A, v, lambda) => {
        const Av = A.dot(Matrix.from([Array.from(v)]).T).col(0);
        return Math.sqrt(Array.from(Av).reduce((s, x, i) => s + (x - lambda * v[i]) ** 2, 0));
    };

    it("should find the smallest eigenpairs of a path graph Laplacian", () => {
        // Like the alignment matrices of LLE and LTSA: positive semidefinite, eigenvalue 0 with the
        // constant vector. The eigenvalues are 2 - 2 cos(jπ / n).
        const n = 60;
        const L = new Matrix(n, n, (i, j) => {
            if (i === j) return i === 0 || i === n - 1 ? 1 : 2;
            return Math.abs(i - j) === 1 ? -1 : 0;
        });
        const { eigenvalues, eigenvectors } = simultaneous_poweriteration(L, 4, { seed: 42, smallest: true });
        for (let j = 0; j < 4; ++j) {
            expect(eigenvalues[j]).toBeCloseTo(2 - 2 * Math.cos((j * Math.PI) / n), 8);
            // the default tol (1e-8 on the squared change per step) leaves the vectors ~1e-5 off
            expect(residual(L, eigenvectors[j], eigenvalues[j])).toBeLessThan(1e-5);
        }
        const c = 1 / Math.sqrt(n);
        for (const x of eigenvectors[0]) expect(Math.abs(x)).toBeCloseTo(c, 6);
    });

    it("should separate tiny, closely spaced eigenvalues", () => {
        // The case power iteration on a shifted matrix cannot handle: the wanted eigenvalues are
        // 1e-9 apart, far below the largest one.
        const n = 40;
        const spectrum = Array.from({ length: n }, (_, i) => (i < 4 ? i * 1e-9 : 1 + i));
        const r = new Randomizer(7);
        const { Q } = qr(new Matrix(n, n, () => r.random - 0.5));
        const A = Q.dot(new Matrix(n, n, (i, j) => (i === j ? spectrum[i] : 0))).dotTrans(Q);
        const { eigenvalues, eigenvectors } = simultaneous_poweriteration(A, 4, { seed: 42, smallest: true });
        for (let j = 0; j < 4; ++j) {
            expect(eigenvalues[j]).toBeCloseTo(spectrum[j], 12);
            expect(residual(A, eigenvectors[j], eigenvalues[j])).toBeLessThan(1e-8);
        }
    });
});
