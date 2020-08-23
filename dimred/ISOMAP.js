import { simultaneous_poweriteration} from "../linear_algebra/index";
import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { Heap } from "../datastructure/index";
import { DR } from "./DR.js";

/**
 * @class
 * @alias ISOMAP
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
        this._k = neighbors || Math.max(Math.floor(X.shape[0] / 10), 2);
    }

    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    transform() {
        let X = this.X;
        let rows = X.shape[0];
        // TODO: make knn extern and parameter for constructor or transform?
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
            kNearestNeighbors.push(H.toArray().slice(1, this._k + 1))
        }
        
        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths
        // TODO: make extern
        /** @see {@link https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm} */
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
        let { eigenvectors: V } = simultaneous_poweriteration(B, this._d);
        this.Y = Matrix.from(V).transpose();
        // return embedding
        return this.projection;
    }


    /**
     * Set and get parameters
     * @param {String} name - name of the parameter.
     * @param {Number} [value = null] - value of the parameter to set, if null then return actual parameter value.
     */
    parameter(name, value=null) {
        return super.parameter(name, value);
    }

    /**
     * Alias for 'parameter'.
     * @param {String} name 
     * @param {Number} value 
     */
    para(name, value=null) {
        return this.parameter(name, value);
    }

    /**
     * Alias for 'parameter'.
     * @param {String} name 
     * @param {Number} value 
     */
    p(name, value=null) {
        return this.parameter(name, value);
    }
}