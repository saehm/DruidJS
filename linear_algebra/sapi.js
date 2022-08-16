import { euclidean } from "../metrics/index.js";
import { Matrix, norm } from "../matrix/index.js";

/**
 * Computes largest Eigenvalue / Eigenvector with
 * Accelerated Stochastic Power Iteration
 * http://proceedings.mlr.press/v84/xu18a/xu18a-supp.pdf
 * @param {Matrix} A - The respective matrix
 * @param {Number} max_iterations - number of maximal iterations
 * @param {Number} batch_size - defines batchsize
 * @param {Number} beta - learning parameter
 * @param {Metric} metric - for computing the norm
 */
export default function (A, max_iterations = 100, batch_size = 10, beta = 0.05, metric = euclidean) {
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    let n = A.shape[0];
    let r = new Matrix(n, 1, () => Math.random());
    let r_last = new Matrix(n, 1, 0);

    while (max_iterations--) {
        let A_r = A.dot(r);
        let beta_r_last = r_last.mult(beta);
        let r_next = A_r.sub(beta_r_last);
        let r_next_norm = norm(r_next._data, metric);
        r_last = r.divide(r_next_norm);
        r = r_next.divide(r_next_norm);
    }

    let u = r.transDot(A).dot(r);
    let l = r.transDot(r);
    let lambda = u.divide(l).entry(0, 0);
    return {
        eigenvector: r.transpose().to2dArray[0],
        eigenvalue: lambda,
    };
}
