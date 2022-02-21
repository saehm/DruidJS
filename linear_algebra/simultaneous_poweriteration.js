import { qr as qr_gramschmidt } from "./index.js";
import { Matrix } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";
import { euclidean_squared } from "../metrics/index.js";

/**
 * Computes the {@link k} biggest Eigenvectors and Eigenvalues from Matrix {@link A} with the QR-Algorithm.
 * @memberof module:linear_algebra
 * @alias simultaneous_poweriteration
 * @param {Matrix} A - The Matrix
 * @param {Number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {Number} [max_iterations=100] - The number of maxiumum iterations the algorithm should run.
 * @param {Number|Randomizer} [seed=1212] - The seed value or a randomizer used in the algorithm.
 * @param {Number} [tol=1e-8] - Allowed error for stopping criteria
 * @returns {{eigenvalues: Array, eigenvectors: Array}} - The {@link k} biggest eigenvectors and eigenvalues of Matrix {@link A}.
 */
export default function (A, k = 2, max_iterations = 100, seed = 1212, qr = qr_gramschmidt, tol = 1e-8) {
    const randomizer = seed instanceof Randomizer ? seed : new Randomizer(seed);
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    const n = A.shape[0];
    let { Q: Q, R: R } = qr(new Matrix(n, k, () => randomizer.random));
    while (max_iterations--) {
        const oldQ = Q.clone();
        const Z = A.dot(Q);
        const QR = qr(Z);
        Q = QR.Q;
        R = QR.R;
        const error = euclidean_squared(Q.values, oldQ.values);
        if (error < tol) {
            break;
        }
    }

    const eigenvalues = R.diag;
    const eigenvectors = Q.transpose().to2dArray;
    return { eigenvalues, eigenvectors };
}
