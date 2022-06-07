import { qr as qr_gramschmidt } from "./index.js";
import { Matrix } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";
import { euclidean_squared } from "../metrics/index.js";

/**
 * Computes the `k` biggest Eigenvectors and Eigenvalues from Matrix `A` with the QR-Algorithm.
 * @memberof module:linear_algebra
 * @alias simultaneous_poweriteration
 * @param {Matrix} A - The Matrix
 * @param {Number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {Object} parameters - Object containing parameterization of the simultanious poweriteration method.
 * @param {Number} [parameters.max_iterations=100] - The number of maxiumum iterations the algorithm should run.
 * @param {Number|Randomizer} [parameters.seed=1212] - The seed value or a randomizer used in the algorithm.
 * @param {Function} [parameters.qr=qr_gramschmidt] - The QR technique to use.
 * @param {Number} [parameters.tol=1e-8] - Tolerated error for stopping criteria.
 * @returns {{eigenvalues: Number[], eigenvectors: Number[][]}} the `k` biggest eigenvectors and eigenvalues of Matrix `A`.
 */
export default function (A, k = 2, {seed = 1212, max_iterations = 100, qr = qr_gramschmidt, tol = 1e-8} = {}) {
    const randomizer = seed instanceof Randomizer ? seed : new Randomizer(seed);
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    const n = A.shape[0];
    let { Q, R } = qr(new Matrix(n, k, () => (randomizer.random - .5) * 2));
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

    const eigenvalues = R.diag;
    const eigenvectors = Q.transpose().to2dArray;
    return { eigenvalues, eigenvectors };
}
