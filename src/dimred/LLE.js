import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { k_nearest_neighbors } from "../matrix/index.js";
import { neumair_sum } from "../numerical/index.js";
import { DR } from "./DR.js";

/**
 * @class
 * @alias LLE
 * @extends DR
 */
export class LLE extends DR {
    /**
     * Locally Linear Embedding.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LLE
     * @param {Matrix} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} parameters.neighbors - the label / class of each data point.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the dimensionality of the projection.
     * @param {object} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://doi.org/10.1126/science.290.5500.2323}
     */
    constructor(X, parameters) {
        super(X, { neighbors: undefined, d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
        this.parameter("neighbors", Math.min(this._parameters.neighbors ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const rows = this._N;
        const cols = this._D;
        const { neighbors, d, eig_args, metric } = this._parameters;
        const nN = k_nearest_neighbors(X, neighbors, metric);
        const O = new Matrix(neighbors, 1, 1);
        const W = new Matrix(rows, rows);

        for (let row = 0; row < rows; ++row) {
            const nN_row = nN[row];
            const Z = new Matrix(neighbors, cols, (i, j) => X.entry(nN_row[i].j, j) - X.entry(row, j));
            const C = Z.dotTrans(Z);
            if (neighbors > cols) {
                const C_trace = neumair_sum(C.diag) / 1000;
                for (let j = 0; j < neighbors; ++j) {
                    C.add_entry(j, j, C_trace);
                }
            }
            // reconstruct;
            let w = Matrix.solve_CG(C, O, this._randomizer);
            w = w.divide(w.sum);
            for (let j = 0; j < neighbors; ++j) {
                W.set_entry(row, nN_row[j].j, w.entry(j, 0));
            }
        }
        // comp embedding
        const I = new Matrix(rows, rows, "identity");
        const IW = I.sub(W);
        const M = IW.transDot(IW);
        const { eigenvectors: V } = simultaneous_poweriteration(M.T.inverse(), d + 1, eig_args);
        this.Y = Matrix.from(V.slice(1, 1 + d)).T;

        // return embedding
        return this.projection;
    }
}
