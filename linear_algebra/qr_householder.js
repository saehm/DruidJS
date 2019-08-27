import { Matrix, norm } from "../matrix/index";

export default function(A) {
    let [rows, cols] = A.shape;
    let Q = new Matrix(rows, rows, "identity");
    let R = A.clone()

    for (let j = 0; j < cols; ++j) {
        let x = R.get_block(j, j).col(0);
        let x_norm = norm(x)
        let x0 = x[0]
        let rho = -Math.sign(x0);
        let u1 = x0 - rho * x_norm;
        let u = x.map(e => e / u1);
        u[0] = 1;
        let beta = -rho * u1 / x_norm;

        //let u_outer_u = new Matrix(cols - j, cols - j, (i, j) => u[i] * u[j]);
        u = Matrix.from(u, "col");
        let R_j_0 = R.get_block(j, 0);
        R.set_block(j, 0, R_j_0.sub(u.dot(u.T.dot(R_j_0)).mult(beta)))
        let Q_0_j = Q.get_block(0, j);
        Q.set_block(0, j, Q_0_j.sub(Q_0_j.dot(u).dot(u.T))) 
    }

    // repair R // numerical unstability?
    for (let i = 1; i < rows; ++i) {
        for (let j = 0; j < i; ++j) {
            R.set_entry(i, j, 0)
        }
    }
    return { R: R, Q: Q };
}

