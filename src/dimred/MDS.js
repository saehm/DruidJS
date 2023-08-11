import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { distance_matrix, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";

/**
 * @class
 * @alias MDS
 * @extends DR
 */
export class MDS extends DR {
    /**
     * Classical MDS.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias MDS
     * @param {Matrix} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {object} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     */
    constructor(X, parameters) {
        super(X, { d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimensionality {@link d}.
     * @returns {Matrix|number[][]}
     */
    transform() {
        const X = this.X;
        const rows = X.shape[0];
        const { d, metric, eig_args } = this._parameters;
        const A = metric === "precomputed" ? X : distance_matrix(X, metric);
        const ai_ = A.meanCols;
        const a_j = A.meanRows;
        const a__ = A.mean;

        this._d_X = A;
        const B = new Matrix(rows, rows, (i, j) => A.entry(i, j) - ai_[i] - a_j[j] + a__);

        const { eigenvectors: V } = simultaneous_poweriteration(B, d, eig_args);
        this.Y = Matrix.from(V).transpose();

        return this.projection;
    }

    /**
     * @returns {number} - the stress of the projection.
     */
    stress() {
        const N = this.X.shape[0];
        const Y = this.Y;
        const d_X = this._d_X;
        const d_Y = new Matrix();
        d_Y.shape = [
            N,
            N,
            (i, j) => {
                return i < j ? euclidean(Y.row(i), Y.row(j)) : d_Y.entry(j, i);
            },
        ];
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
