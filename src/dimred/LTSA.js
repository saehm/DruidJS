import { Matrix, k_nearest_neighbors } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { DR } from "./DR.js";

/**
 * @class
 * @alias LTSA
 * @extends DR
 */
export class LTSA extends DR {
    /**
     * Local Tangent Space Alignment
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LTSA
     * @param {Matrix} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} parameters.neighbors - the number of neighbors {@link LTSA} should use to project the data.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {object} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
     */
    constructor(X, parameters) {
        super(X, { neighbors: undefined, d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
        this.parameter("neighbors", Math.min(this._parameters.neighbors ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
        if (this._D <= this.parameter("d")) {
            throw new Error(`Dimensionality of X (D = ${this._D}) must be greater than the required dimensionality of the result (d = ${this.parameter("d")})!`);
        }
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const [rows, D] = X.shape;
        const { d, neighbors, metric, eig_args } = this._parameters;
        // 1.1 determine k nearest neighbors
        const nN = k_nearest_neighbors(X, neighbors, metric);
        // center matrix
        const O = new Matrix(D, D, "center");
        const B = new Matrix(rows, rows, 0);

        for (let row = 0; row < rows; ++row) {
            // 1.2 compute the d largest eigenvectors of the correlation matrix
            const I_i = [row, ...nN[row].map((n) => n.j)];
            let X_i = Matrix.from(I_i.map((n) => X.row(n)));
            // center X_i
            X_i = X_i.dot(O);
            // correlation matrix
            const C = X_i.dotTrans(X_i);
            const { eigenvectors: g } = simultaneous_poweriteration(C, d, eig_args);
            //g.push(linspace(0, k).map(_ => 1 / Math.sqrt(k + 1)));
            const G_i_t = Matrix.from(g);
            // 2. Constructing alignment matrix
            const W_i = G_i_t.transDot(G_i_t).add(1 / Math.sqrt(neighbors + 1));
            for (let i = 0; i < neighbors + 1; ++i) {
                for (let j = 0; j < neighbors + 1; ++j) {
                    B.add_entry(I_i[i], I_i[j], W_i.entry(i, j) - (i === j ? 1 : 0));
                }
            }
        }

        // 3. Aligning global coordinates
        const { eigenvectors: Y } = simultaneous_poweriteration(B, d + 1, eig_args);
        this.Y = Matrix.from(Y.slice(1)).transpose();

        // return embedding
        return this.projection;
    }
}
