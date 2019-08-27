import { simultaneous_poweriteration} from "../linear_algebra/index"
//import { dijkstra } from "../knn/index";
import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { Heap } from "../datastructure/index";

export class ISOMAP{
    constructor(X, neighbors, d=2, metric=euclidean) {
        this.X = X;
        this.k = neighbors || Math.floor(this.X.shape[0] / 10)
        this.d = d;
        this._metric = metric;
    }

    transform() {
        let X = this.X;
        let rows = X.shape[0];
        // make knn extern and parameter for constructor or transform?
        let D = new Matrix()
        D.shape = [rows, rows, (i,j) => i <= j ? this._metric(X.row(i), X.row(j)) : D.entry(j,i)]
        let kNearestNeighbors = [];
        for (let i = 0; i < rows; ++i) {
            let row = D.row(i).map((d,i) => { 
                return {
                    "index": i,
                    "distance": d
                }
            });
            let H = new Heap(row, d => d.distance, "min");
            kNearestNeighbors.push(H.toArray().slice(1, this.k + 1))
        }
        
        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths
        // TODO: make extern
        let G = new Matrix(rows, rows, (i,j) => {
            let other = kNearestNeighbors[i].find(n => n.index === j);
            return other ? other.distance : Infinity
        });

        for (let i = 0; i < rows; ++i) {
            for (let j = 0; j < rows; ++j) {
                for (let k = 0; k < rows; ++k) {
                    G.set_entry(i, j, Math.min(G.entry(i, j), G.entry(i, k) + G.entry(k, j)));
                }
            }
        }
        
        let ai_ = [];
        let a_j = [];
        for (let i = 0; i < rows; ++i) {
            ai_.push(0)
            a_j.push(0)
        }
        let a__ = 0;
        let A = new Matrix(rows, rows, (i,j) => {
            let val = G.entry(i, j);
            val = val === Infinity ? 0 : val;
            ai_[i] += val;
            a_j[j] += val;
            a__ += val;
            return val;
        });
        
        ai_ = ai_.map(v => v / rows);
        a_j = a_j.map(v => v / rows);
        a__ /= (rows ** 2);
        let B = new Matrix(rows, rows, (i,j) => (A.entry(i,j) - ai_[i] - a_j[j] + a__));
             
        // compute d eigenvectors
        let { eigenvectors: V } = simultaneous_poweriteration(B, this.d);
        this.Y = Matrix.from(V).transpose();
        // return embedding
        return this.Y
    }

    get projection() {
        return this.Y
    }

    

}