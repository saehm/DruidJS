import { qr } from "./index";
import { Matrix } from "../matrix/Matrix";
import { Randomizer } from "../util/index";
import { neumair_sum } from "../numerical/index";

/**
 * Computes the {@link k} biggest Eigenvectors and Eigenvalues from Matrix {@link A} with the QR-Algorithm.
 * @param {Matrix} A - The Matrix
 * @param {Number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {Number} [max_iterations=100] - The number of maxiumum iterations the algorithm should run.
 * @param {Number|Randomizer} [seed=1212] - The seed value or a randomizer used in the algorithm.
 * @returns {{eigenvalues: Array, eigenvectors: Array}} - The {@link k} biggest eigenvectors and eigenvalues of Matrix {@link A}.
 */
export default function(A, k = 2, max_iterations=100, seed=1212) {
    const randomizer = seed instanceof Randomizer ? seed : new Randomizer(seed);
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    const n = A.shape[0]
    let { Q: Q, R: R } = qr(new Matrix(n, k, () => randomizer.random));
    while (max_iterations--) {
        const oldR = R.clone();
        const Z = A.dot(Q);
        const QR = qr(Z); 
        Q = QR.Q;
        R = QR.R;
        if (neumair_sum(R.sub(oldR).diag) / n < 1e-12) {
            max_iterations = 0;
        }        
    }

    const eigenvalues = R.diag;
    const eigenvectors = Q.transpose().to2dArray;
    return {
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
    };
}

