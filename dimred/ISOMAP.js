import { simultaneous_poweriteration} from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { Heap } from "../datastructure/index.js";
import { DR } from "./DR.js";

/**
 * @class
 * @alias ISOMAP
 * @extends DR
 */
export class ISOMAP extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias ISOMAP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} neighbors - the number of neighbors {@link ISOMAP} should use to project the data.
     * @param {Number} [d = 2] - the dimensionality of the projection. 
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points. 
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
    constructor(X, neighbors, d = 2, metric = euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["k"];
        this.parameter("k", Math.min(neighbors ?? Math.max(Math.floor(this.X.shape[0] / 10), 2), this._N -1));
        return this;
    }

    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    transform() {
        this.check_init();
        const X = this.X;
        const rows = this._N;
        const metric = this._metric;
        // TODO: make knn extern and parameter for constructor or transform?
        const D = new Matrix();
        D.shape = [rows, rows, (i,j) => i <= j ? metric(X.row(i), X.row(j)) : D.entry(j,i)]
        const kNearestNeighbors = [];
        for (let i = 0; i < rows; ++i) {
            const row = [];
            for (let j = 0; j < rows; ++j) {
                row.push({
                    "index": j,
                    "distance": D.entry(i, j),
                })
            }
            const H = new Heap(row, d => d.distance, "min");
            kNearestNeighbors.push(H.toArray().slice(1, this._k + 1))
        }
        
        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths
        // TODO: make extern
        /** @see {@link https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm} */
        const G = new Matrix(rows, rows, (i,j) => {
            const other = kNearestNeighbors[i].find(n => n.index === j);
            return other ? other.distance : Infinity
        });

        for (let i = 0; i < rows; ++i) {
            for (let j = 0; j < rows; ++j) {
                for (let k = 0; k < rows; ++k) {
                    G.set_entry(i, j, Math.min(G.entry(i, j), G.entry(i, k) + G.entry(k, j)));
                }
            }
        }
        
        let ai_ = new Float64Array(rows);
        let a_j = new Float64Array(rows);
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
        const B = new Matrix(rows, rows, (i,j) => (A.entry(i,j) - ai_[i] - a_j[j] + a__));
             
        // compute d eigenvectors
        const { eigenvectors: V } = simultaneous_poweriteration(B, this._d);
        this.Y = Matrix.from(V).transpose();
        // return embedding
        return this.projection;
    }


}