/**
 * Computes the {@link k} biggest Eigenvectors and Eigenvalues from Matrix {@link A} with the QR-Algorithm.
 * @memberof module:linear_algebra
 * @alias simultaneous_poweriteration
 * @param {Matrix} A - The Matrix
 * @param {Number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {Object} parameters - Object containing parameterization of the simultanious poweriteration method.
 * @param {Number} [parameters.max_iterations=100] - The number of maxiumum iterations the algorithm should run.
 * @param {Number|Randomizer} [parameters.seed=1212] - The seed value or a randomizer used in the algorithm.
 * @param {Function} [parameters.qr=qr_gramschmidt] - The QR technique to use.
 * @param {Number} [parameters.tol=1e-8] - Allowed error for stopping criteria
 * @returns {{eigenvalues: Array, eigenvectors: Array}} - The {@link k} biggest eigenvectors and eigenvalues of Matrix {@link A}.
 */
export default function _default(A: Matrix, k?: number, { seed, max_iterations, qr, tol }?: {
    max_iterations?: number;
    seed?: number | Randomizer;
    qr?: Function;
    tol?: number;
}): {
    eigenvalues: any[];
    eigenvectors: any[];
};
import { Matrix } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";
//# sourceMappingURL=simultaneous_poweriteration.d.ts.map