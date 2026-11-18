import { describe, expect, it } from "vitest";
import { inner_product, qr, qr_householder, simultaneous_poweriteration } from "../../src/linear_algebra/index.js";
import { Matrix } from "../../src/matrix/index.js";

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
