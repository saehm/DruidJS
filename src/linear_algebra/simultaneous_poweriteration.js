import { Matrix } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";
import { qr as qr_gramschmidt } from "./index.js";

/** @import { EigenArgs } from "./index.js" */

/**
 * Computes the `k` biggest Eigenvectors and Eigenvalues from Matrix `A` with the QR-Algorithm.
 *
 * "Biggest" means largest in magnitude. The eigenvalues keep their sign.
 *
 * @category Linear Algebra
 * @param {Matrix} A - The Matrix
 * @param {number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {EigenArgs} parameters - Object containing parameterization of the simultanious
 *   poweriteration method.
 * @returns {{ eigenvalues: Float64Array; eigenvectors: Float64Array[] }} The `k` biggest eigenvectors and eigenvalues
 *   of Matrix `A`.
 */
export function simultaneous_poweriteration(
    A,
    k = 2,
    { seed = 1212, max_iterations = 100, qr = qr_gramschmidt, tol = 1e-8 } = {},
) {
    const randomizer = seed instanceof Randomizer ? seed : new Randomizer(seed);
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    const n = A.shape[0];
    let { Q } = qr(new Matrix(n, k, () => (randomizer.random - 0.5) * 2));
    while (max_iterations--) {
        const oldQ = Q;
        const Z = A.dot(Q);
        Q = qr(Z).Q;
        if (converged(Q, oldQ, tol)) {
            break;
        }
    }

    // Rayleigh quotients qᵀAq. Unlike the diagonal of R, which holds only |λ|, they keep the sign
    // of negative eigenvalues (e.g. of the double-centered matrix of non-Euclidean distances).
    const AQ = A.dot(Q);
    const eigenvalues = new Float64Array(k);
    for (let j = 0; j < k; ++j) {
        let lambda = 0;
        for (let i = 0; i < n; ++i) lambda += Q.entry(i, j) * AQ.entry(i, j);
        eigenvalues[j] = lambda;
    }
    const eigenvectors = Q.transpose().to2dArray();
    return { eigenvalues, eigenvectors };
}

/**
 * Whether the columns of `Q` equal those of `oldQ` up to sign. A column whose eigenvalue is negative
 * flips its sign in every step, so comparing the columns directly would never report convergence.
 *
 * @private
 * @param {Matrix} Q
 * @param {Matrix} oldQ
 * @param {number} tol - Bound on the summed squared differences.
 * @returns {boolean}
 */
function converged(Q, oldQ, tol) {
    const [n, k] = Q.shape;
    const q = Q.values;
    const o = oldQ.values;
    let error = 0;
    for (let j = 0; j < k; ++j) {
        let same = 0;
        let flipped = 0;
        for (let i = 0; i < n; ++i) {
            const a = q[i * k + j];
            const b = o[i * k + j];
            same += (a - b) ** 2;
            flipped += (a + b) ** 2;
        }
        error += Math.min(same, flipped);
    }
    return error < tol;
}
