import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { Randomizer } from "../util/index";
import { simultaneous_poweriteration} from "../linear_algebra/index";
import { k_nearest_neighbors } from "../matrix/index";
import { neumair_sum } from "../numerical/index";

export class LLE{
    constructor(X, neighbors, d=2, metric=euclidean) {
        this.X = X;
        this._k = neighbors;
        this.d = d;
        this._metric = metric;
    }

    transform() {
        let X = this.X;
        let d = this.d;
        let [ rows, cols ] = X.shape;
        let k = this._k;
        let nN = k_nearest_neighbors(X.to2dArray, k, null, this._metric);
        let O = new Matrix(k, 1, 1);
        let W = new Matrix(rows, rows);

        for (let row = 0; row < rows; ++row) {
            let Z = new Matrix(k, cols, (i, j) => X.entry(nN[row][i].j, j) - X.entry(row, j));
            let C = Z.dot(Z.transpose());
            if ( k > cols ) {
                let C_trace = neumair_sum(C.diag) / 1000;
                for (let j = 0; j < k; ++j) {
                    C.set_entry(j, j, C.entry(j, j) + C_trace);
                }
            }

            // reconstruct;
            let w = Matrix.solve(C, O);
            let w_sum = neumair_sum(w.col(0));
            w = w.divide(w_sum);
            for (let j = 0; j < k; ++j) {
                W.set_entry(row, nN[row][j].j, w.entry(j, 0));
            }
        }
        // comp embedding
        let I = new Matrix(rows, rows, "identity");
        let IW = I.sub(W);
        let M = IW.transpose().dot(IW);
        let { eigenvectors: V } = simultaneous_poweriteration(M.transpose().inverse(), d + 1);
        
        this.Y = Matrix.from(V.slice(1, 1 + d)).transpose();

        // return embedding
        return this.Y;
    }

    get projection() {
        return this.Y;
    }
}