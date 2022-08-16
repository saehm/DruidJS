import { Matrix, norm } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { neumair_sum } from "../numerical/index.js";

/**
 * Computes the QR Decomposition of the Matrix `A` using Gram-Schmidt process.
 * @memberof module:linear_algebra
 * @alias qr
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process}
 */
export default function (A) {
    const [rows, cols] = A.shape;
    const Q = new Matrix(rows, cols, "identity");
    const R = new Matrix(cols, cols, 0);

    for (let j = 0; j < cols; ++j) {
        let v = A.col(j);
        for (let i = 0; i < j; ++i) {
            const q = Q.col(i);
            const q_dot_v = neumair_sum(q.map((q_, k) => q_ * v[k]));
            for (let k = 0; k < rows; ++k) {
                v[k] -= q_dot_v * q[k];
            }
            R.set_entry(i, j, q_dot_v);
        }
        const v_norm = norm(v, euclidean);
        for (let k = 0; k < rows; ++k) {
            Q.set_entry(k, j, v[k] / v_norm);
        }
        R.set_entry(j, j, v_norm);
    }
    return { R, Q };
}
