import { Matrix, k_nearest_neighbors } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { simultaneous_poweriteration} from "../linear_algebra/index";
import { DR } from "./DR.js";

/**
 * @class
 * @alias LTSA
 */
export class LTSA extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LTSA
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
     */
    constructor(X, neighbors, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["k"];
        this.parameter("k", Math.min(neighbors ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        if (this._D <= d) throw `Dimensionality of X (D = ${this._D}) must be greater than the required dimensionality of the result (d = ${d})!`;
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const d = this._d;
        const [ rows, D ] = X.shape;
        const k = this._k;
        // 1.1 determine k nearest neighbors
        const nN = k_nearest_neighbors(X.to2dArray, k, null, this._metric);
        // center matrix
        const O = new Matrix(D, D, "center");
        const B = new Matrix(rows, rows, 0);
        
        for (let row = 0; row < rows; ++row) {
            // 1.2 compute the d largest eigenvectors of the correlation matrix
            const I_i = [row, ...nN[row].map(n => n.j)]
            let X_i = Matrix.from(I_i.map(n => X.row(n)));
            // center X_i
            X_i = X_i.dot(O)
            // correlation matrix
            const C = X_i.dot(X_i.transpose());
            const { eigenvectors: g } = simultaneous_poweriteration(C, d);
            //g.push(linspace(0, k).map(_ => 1 / Math.sqrt(k + 1)));
            const G_i_t = Matrix.from(g);
            // 2. Constructing alignment matrix
            const W_i = G_i_t.transpose().dot(G_i_t).add(1 / Math.sqrt(k + 1));
            for (let i = 0; i < k + 1; ++i) {
                for (let j = 0; j < k + 1; ++j) {
                    B.set_entry(I_i[i], I_i[j], B.entry(I_i[i], I_i[j]) - (i === j ? 1 : 0 ) + W_i.entry(i, j));
                }
            }
        }

        // 3. Aligning global coordinates
        const { eigenvectors: Y } = simultaneous_poweriteration(B, d + 1);
        this.Y = Matrix.from(Y.slice(1)).transpose();

        // return embedding
        return this.projection;
    }
}