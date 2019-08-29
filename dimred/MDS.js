import { simultaneous_poweriteration} from "../linear_algebra/index";
import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";

export class MDS{
    constructor(X, d=2, metric=euclidean) {
        this.X = X;
        this.d = d;
        this._metric = metric;
    }

    transform() {
        const X = this.X;
        //let sum_reduce = (a,b) => a + b
        const rows = X.shape[0];
        const metric = this._metric;
        let ai_ = [];
        let a_j = [];
        for (let i = 0; i < rows; ++i) {
            ai_.push(0)
            a_j.push(0)
        }
        let a__ = 0;
        const A = new Matrix();
        A.shape = [rows, rows, (i,j) => {
            let val = 0
            if (i < j) {
                val = metric(X.row(i), X.row(j))
            } else if (i > j) {
                val = A.entry(j,i);
            }
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
                
        const { eigenvectors: V } = simultaneous_poweriteration(B, this.d);
        this.Y = Matrix.from(V).transpose()
        
        return this.Y
    }

    get projection() {
        return this.Y
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