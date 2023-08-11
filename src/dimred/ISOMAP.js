import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";
import { Heap } from "../datastructure/index.js";
import { DR } from "./DR.js";
import euclidean from "../metrics/euclidean.js";

/**
 * @class
 * @alias ISOMAP
 * @extends DR
 */
export class ISOMAP extends DR {
    /**
     * Isometric feature mapping (ISOMAP).
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias ISOMAP
     * @param {Matrix} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} parameters.neighbors - the number of neighbors {@link ISOMAP} should use to project the data.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {object} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://doi.org/10.1126/science.290.5500.2319}
     */
    constructor(X, parameters) {
        super(X, { neighbors: undefined, d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
        this.parameter("neighbors", Math.min(this._parameters.neighbors ?? Math.max(Math.floor(this.X.shape[0] / 10), 2), this._N - 1));
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
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
        const { d, metric, eig_args, neighbors } = this._parameters;
        // TODO: make knn extern and parameter for constructor or transform?
        const D = new Matrix();
        D.shape = [rows, rows, (i, j) => (i <= j ? metric(X.row(i), X.row(j)) : D.entry(j, i))];
        const kNearestNeighbors = [];
        for (let i = 0; i < rows; ++i) {
            const row = [];
            for (let j = 0; j < rows; ++j) {
                row.push({
                    index: j,
                    distance: D.entry(i, j),
                });
            }
            const H = new Heap(row, (d) => d.distance, "min");
            kNearestNeighbors.push(H.toArray().slice(1, neighbors + 1));
        }

        /*D = dijkstra(kNearestNeighbors);*/
        // compute shortest paths
        // TODO: make extern
        /** @see {@link https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm} */
        const G = new Matrix(rows, rows, (i, j) => {
            const other = kNearestNeighbors[i].find((n) => n.index === j);
            return other ? other.distance : Infinity;
        });

        for (let i = 0; i < rows; ++i) {
            for (let j = 0; j < rows; ++j) {
                let min_val = G.entry(i, j);
                for (let k = 0; k < rows; ++k) {
                    min_val = Math.min(min_val, G.entry(i, k) + G.entry(k, j));
                }
                G.set_entry(i, j, min_val);
            }
        }

        let ai_ = new Float64Array(rows);
        let a_j = new Float64Array(rows);
        let a__ = 0;
        const A = new Matrix(rows, rows, (i, j) => {
            let val = G.entry(i, j);
            val = val === Infinity ? 0 : val;
            ai_[i] += val;
            a_j[j] += val;
            a__ += val;
            return val;
        });

        ai_ = ai_.map((v) => v / rows);
        a_j = a_j.map((v) => v / rows);
        a__ /= rows ** 2;
        const B = new Matrix(rows, rows, (i, j) => A.entry(i, j) - ai_[i] - a_j[j] + a__);

        // compute d eigenvectors
        const { eigenvectors: V } = simultaneous_poweriteration(B, d, eig_args);
        this.Y = Matrix.from(V).transpose();
        // return embedding
        return this.projection;
    }
}
