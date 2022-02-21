import { Matrix, norm } from "../matrix/index.js";

/**
 * Computes the QR Decomposition of the Matrix {@link A} with householder transformations.
 * @memberof module:linear_algebra
 * @alias qr_householder
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections}
 * @see {@link http://mlwiki.org/index.php/Householder_Transformation}
 */
export default function (A) {
    const [rows, cols] = A.shape;
    const Q = new Matrix(rows, rows, "I");
    const R = A.clone();

    for (let j = 0; j < cols; ++j) {
        const x = Matrix.from(R.col(j).slice(j));
        const x_norm = norm(x);
        const x0 = x.entry(0, 0);
        const rho = -Math.sign(x0);
        const u1 = x0 - rho * x_norm;
        const u = x.divide(u1).set_entry(0, 0, 1);
        const beta = (-rho * u1) / x_norm;

        const u_outer_u = u.outer(u);
        const R_block = R.get_block(j, 0);
        const new_R = R_block.sub(u_outer_u.dot(R_block).mult(beta));
        const Q_block = Q.get_block(0, j);
        const new_Q = Q_block.sub(Q_block.dot(u_outer_u).mult(beta));
        R.set_block(j, 0, new_R);
        Q.set_block(0, j, new_Q);
    }
    return { R, Q };
}
