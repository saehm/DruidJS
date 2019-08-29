import { Matrix, norm } from "../matrix/index";
import { euclidean } from "../metrics/index"
import { neumair_sum } from "../numerical/index";

/**
 * Computes the QR Decomposition of the Matrix {@link A}.
 * @memberof module:linear_algebra
 * @alias qr
 * @param {Matrix} A 
 */
export default function(A) {
    let [rows, cols] = A.shape;
    let Q = new Matrix(rows, cols, "identity");
    let R = new Matrix(cols, cols, 0);

    for (let j = 0; j < cols; ++j) {
        let v = A.col(j);
        for (let i = 0; i < j; ++i) {
            let q = Q.col(i);
            let q_dot_v = neumair_sum(q.map((q_, k) => q_ * v[k]));
            R.set_entry(i,j, q_dot_v);
            v = v.map((v_, k) => v_ - q_dot_v * q[k]);
        }
        let v_norm = norm(v, euclidean);
        for (let k = 0; k < rows; ++k) {
            Q.set_entry(k, j, v[k] / v_norm);
        }
        R.set_entry(j,j, v_norm)
    }
    /*let R = matrix.clone()

    for (let j = 0; j < cols; ++j) {
        let z = new Matrix()
        z.shape = [rows - j, 1, (i,_) => R.entry(i+j, j)]
        let norm_z = norm(z, euclidean);
        let rho = z.entry(0, 0) < 0 ? 1 : -1;
        let u1 = z.entry(0, 0) - rho * norm_z;
        let u = z.divide(u1);
        u.set_entry(0, 0, 1);
        let beta = -rho * u1 / norm_z;

        let row_j = rows - j
        let u_outer_u = u.dot(u.transpose())
        let R_j__ = new Matrix()
        R_j__.shape = [row_j, row_j, (r, c) => R.entry(row_j + r, row_j + c)]
        let R_sub = u_outer_u.dot(R_j__).mult(-beta)
        for (let r = 0; r < row_j; ++r) {
            for (let c = 0; c < row_j; ++c) {
                R.set_entry(row_j + r, row_j + c, R.entry(row_j + r, row_j + c) + R_sub.entry(r, c))
            }
        }
        let Q__j_ = new Matrix();
        Q__j_.shape = [rows, row_j, (r, c) => Q.entry(r, row_j + c)]
        let Q_sub = Q__j_.dot(u_outer_u).mult(-beta)
        for (let r = 0; r < rows; ++r) {
            for ( let c = 0; c < row_j; ++c) {
                Q.set_entry(r, row_j + c, Q.entry(r, row_j + c) + Q_sub.entry(r, c))
            }
        }

        console.log(j, Q, R)
    }*/

    return { R: R, Q: Q };
}

