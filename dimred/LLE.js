import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { simultaneous_poweriteration} from "../linear_algebra/index";
import { k_nearest_neighbors } from "../matrix/index";
import { neumair_sum } from "../numerical/index";
import { DR } from "./DR.js";

/**
 * @class
 * @alias LLE
 */
export class LLE extends DR {
    static parameter_list = ["k"];

    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LLE
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
    constructor(X, neighbors, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = druid.LLE.parameter_list;
        this.parameter("k", neighbors);
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const d = this._d;
        const [ rows, cols ] = X.shape;
        const k = this._k;
        const nN = k_nearest_neighbors(X.to2dArray, k, null, this._metric);
        const O = new Matrix(k, 1, 1);
        const W = new Matrix(rows, rows);

        for (let row = 0; row < rows; ++row) {
            const Z = new Matrix(k, cols, (i, j) => X.entry(nN[row][i].j, j) - X.entry(row, j));
            const C = Z.dot(Z.transpose());
            if ( k > cols ) {
                const C_trace = neumair_sum(C.diag) / 1000;
                for (let j = 0; j < k; ++j) {
                    C.set_entry(j, j, C.entry(j, j) + C_trace);
                }
            }

            // reconstruct;
            let w = Matrix.solve(C, O);
            const w_sum = neumair_sum(w.col(0));
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
        return this.projection;
    }
}