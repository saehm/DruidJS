import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { k_nearest_neighbors } from "../matrix/index.js";
import { neumair_sum } from "../numerical/index.js";
import { norm } from "../matrix/index.js";

export class MLLE {
    constructor(X, neighbors, d = 2, metric = euclidean) {
        this.X = X;
        this._k = neighbors;
        this.d = d;
        this._metric = metric;
    }

    transform() {
        let X = this.X;
        let d = this.d;
        let [rows, cols] = X.shape;
        let k = this._k;
        // 1.1 Determine a neighborset
        let nN = k_nearest_neighbors(X.to2dArray, k, null, this._metric);
        let O = new Matrix(k, 1, 1);
        let W = new Matrix(rows, k);
        let Phi = new Matrix(rows, rows);

        let V = new Array(rows);
        let Lambda = new Array(rows);
        let P = new Array(rows);

        for (let row = 0; row < rows; ++row) {
            let I_i = nN[row].map((n) => n.j);
            let x_i = Matrix.from(X.row(row), "row");
            let X_i = Matrix.from(I_i.map((n) => X.row(n)));
            X_i = X_i.sub(x_i);
            //X_i = X_i.dot(new Matrix(X_i._cols, X_i._cols, "center"))
            let C_i = X_i.dotTrans(X_i); // k by k

            let gamma = neumair_sum(C_i.diag) / 1000;
            for (let j = 0; j < k; ++j) {
                C_i.set_entry(j, j, C_i.entry(j, j) + gamma);
            }

            let { eigenvalues: Lambdas, eigenvectors: v } = simultaneous_poweriteration(C_i, k);
            V[row] = v; // k by k, rows are eigenvectors, big to small
            Lambda[row] = Lambdas; // 1 by k, cols are eigenvalues, big to small
            P.push(neumair_sum(Lambdas.slice(d + 1)) / neumair_sum(Lambdas.slice(0, d)));

            // reconstruct;
            let w = Matrix.solve(C_i, O); // k by 1
            let w_sum = neumair_sum(w.col(0));
            w = w.divide(w_sum);
            for (let j = 0; j < k; ++j) {
                W.set_entry(row, j, w.entry(j, 0));
            }
        }
        // find regularized weights // median
        let theta = P.sort((rho_i, rho_j) => rho_i - rho_j)[Math.ceil(rows / 2)];

        for (let row = 0; row < rows; ++row) {
            let I_i = nN[row].map((n) => n.j);
            let Lambdas = Lambda[row]; // 1 by k
            let s_i = Lambdas.map((Lambda, l) => {
                return {
                    l: l,
                    ratio: neumair_sum(Lambdas.slice(k - l + 1)) / neumair_sum(Lambdas.slice(0, k - l)),
                    Lambda: Lambda,
                };
            });
            //console.log(s_i)
            s_i =
                s_i
                    .filter((s) => s.ratio < theta && s.l <= k - d)
                    .map((s) => s.l)
                    .pop() || d;
            let V_i = V[row]; // k by k
            V_i = V_i.slice(k - s_i); // s_i by k
            let alpha_i = (1 / Math.sqrt(s_i)) * norm(V_i[0].map((_, j) => neumair_sum(V_i.map((r) => r[j]))));
            V_i = Matrix.from(V_i); // s_i by k

            //https://github.com/scikit-learn/scikit-learn/blob/7b136e9/sklearn/manifold/locally_linear.py#L703

            let h = new Matrix(s_i, 1, alpha_i);
            let ones = new Matrix(k, 1, 1);
            h = h.sub(V_i.dot(ones));
            let h_norm = norm(h.col(0));
            h = h_norm < 1e-12 ? h.mult(0) : h.divide(h_norm);
            V_i = V_i.T;
            ones = new Matrix(s_i, 1, 1);
            let w_i = Matrix.from(W.row(row), "col");

            /*let H_i = new Matrix(s_i, s_i, "identity");
            H_i = H_i.sub(h.mult(2).outer(h));
            let W_i = V_i.sub(V_i.dot(h).dot(h.T).mult(2)).add(w_i.mult(1 - alpha_i))
            */
            let W_i = V_i.sub(V_i.dot(h).dotTrans(h).mult(2)).add(w_i.mult(1 - alpha_i).dotTrans(ones));

            W_i = W_i.dotTrans(W_i);
            for (let i = 0; i < k + 1; ++i) {
                for (let j = 0; j < s_i; ++j) {
                    Phi.set_entry(I_i[i], I_i[j], Phi.entry(I_i[i], I_i[j]) - (i === j ? 1 : 0) + W_i.entry(i, j));
                }
            }
        }
        //let { eigenvectors: Y } = simultaneous_poweriteration(Phi.inverse(), d + 1);
        //this.Y = Matrix.from(Y.slice(1)).transpose()

        let { eigenvectors: Y } = simultaneous_poweriteration(Phi, d + 1);
        this.Y = Matrix.from(Y.slice(1)).transpose();

        // return embedding
        return this.Y;
    }

    get projection() {
        return this.Y;
    }
}
