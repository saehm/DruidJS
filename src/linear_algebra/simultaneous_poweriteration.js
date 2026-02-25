import { Matrix } from "../matrix/index.js";
import { euclidean_squared } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { qr as qr_gramschmidt } from "./index.js";

/** @import { EigenArgs } from "./index.js" */

/**
 * Computes the `k` biggest Eigenvectors and Eigenvalues from Matrix `A` with the QR-Algorithm.
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
    let { Q, R } = qr(new Matrix(n, k, () => (randomizer.random - 0.5) * 2));
    while (max_iterations--) {
        const oldQ = Q;
        const Z = A.dot(Q);
        const QR = qr(Z);
        Q = QR.Q;
        R = QR.R;
        const error = euclidean_squared(Q.values, oldQ.values);
        if (error < tol) {
            break;
        }
    }

    const eigenvalues = R.diag();
    const eigenvectors = Q.transpose().to2dArray();
    return { eigenvalues, eigenvectors };
}
