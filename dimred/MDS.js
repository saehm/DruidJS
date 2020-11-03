import { simultaneous_poweriteration} from "../linear_algebra/index";
import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { DR } from "./DR.js";

/**
 * @class
 * @alias MDS
 */
export class MDS extends DR{
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias MDS
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
    
    constructor(X, d=2, metric=euclidean, seed=1212) {
        super(X, d, metric, seed);
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const rows = X.shape[0];
        const metric = this._metric;
        let ai_ = new Float64Array(rows);
        let a_j = new Float64Array(rows);
        let a__ = 0;

        const A = new Matrix();
        A.shape = [rows, rows, (i,j) => {
            if (i === j) return 0;
            const val = (i < j) ? metric(X.row(i), X.row(j)) : A.entry(j,i);
            ai_[i] += val;
            a_j[j] += val;
            a__ += val;
            return val;
        }];
        this._d_X = A;
        ai_ = ai_.map(v => v / rows);
        a_j = a_j.map(v => v / rows);
        a__ /= (rows ** 2);
        const B = new Matrix(rows, rows, (i, j) => (A.entry(i, j) - ai_[i] - a_j[j] + a__));
        //B.shape = [rows, rows, (i,j) => (A.entry(i,j) - (A.row(i).reduce(sum_reduce) / rows) - (A.col(j).reduce(sum_reduce) / rows) + a__)]
                
        const { eigenvectors: V } = simultaneous_poweriteration(B, this._d);
        this.Y = Matrix.from(V).transpose()
        
        return this.projection;
    }

    get stress() {
        const N = this.X.shape[0];
        const Y = this.Y;
        const d_X = this._d_X; /*new Matrix();
        d_X.shape = [N, N, (i, j) => {
            return i < j ? metric(X.row(i), X.row(j)) : d_X.entry(j, i);
        }]*/
        const d_Y = new Matrix();
        d_Y.shape = [N, N, (i, j) => {
            return i < j ? euclidean(Y.row(i), Y.row(j)) : d_Y.entry(j, i);
        }]
        let top_sum = 0;
        let bottom_sum = 0;
        for (let i = 0; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                top_sum += Math.pow(d_X.entry(i, j) - d_Y.entry(i, j), 2);
                bottom_sum += Math.pow(d_X.entry(i, j), 2);
            }
        }
        return Math.sqrt(top_sum / bottom_sum);
    }
}