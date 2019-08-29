import { Matrix, k_nearest_neighbors, linspace } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { Randomizer } from "../util/randomizer";
import { simultaneous_poweriteration} from "../linear_algebra/index";

// https://epubs.siam.org/doi/abs/10.1137/S1064827502419154
export class LTSA{

    constructor(X, neighbors, d=2, metric=euclidean) {
        this.X = X;
        this._k = neighbors;
        this.d = d;
        this._metric = metric;
    }

    transform() {
        let X = this.X;
        let d = this.d;
        let [ rows, D ] = X.shape;
        let k = this._k;
        // 1.1 determine k nearest neighbors
        let nN = k_nearest_neighbors(X.to2dArray, k, null, this._metric);
        // center matrix
        let O = new Matrix(D, D, "center");
        let B = new Matrix(rows, rows, 0);
        
        for (let row = 0; row < rows; ++row) {
            // 1.2 compute the d largest eigenvectors of the correlation matrix
            let I_i = [row, ...nN[row].map(n => n.j)]
            let X_i = Matrix.from(I_i.map(n => X.row(n)));
            // center X_i
            X_i = X_i.dot(O)
            // correlation matrix
            let C = X_i.dot(X_i.transpose());
            let { eigenvectors: g } = simultaneous_poweriteration(C, d);
            //g.push(linspace(0, k).map(_ => 1 / Math.sqrt(k + 1)));
            let G_i_t = Matrix.from(g);
            // 2. Constructing alignment matrix
            let W_i = G_i_t.transpose().dot(G_i_t).add(1 / Math.sqrt(k + 1));
            for (let i = 0; i < k + 1; ++i) {
                for (let j = 0; j < k + 1; ++j) {
                    B.set_entry(I_i[i], I_i[j], B.entry(I_i[i], I_i[j]) - (i === j ? 1 : 0 ) + W_i.entry(i, j));
                }
            }
        }

        // 3. Aligning global coordinates
        let { eigenvectors: Y } = simultaneous_poweriteration(B, d + 1);
        this.Y = Matrix.from(Y.slice(1)).transpose()

        // return embedding
        return this.Y;
    }

    get projection() {
        return this.Y;
    }
}