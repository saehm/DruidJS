import { spatial_tree } from "../knn/index.js";
import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersLTSA} from "./index.js" */
/** @import {KNN} from "../knn/KNN.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */

/**
 * Local Tangent Space Alignment (LTSA)
 *
 * A nonlinear dimensionality reduction algorithm that represents the local
 * geometry of the manifold by tangent spaces and then aligns them to reveal
 * the global structure.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersLTSA>
 * @category Dimensionality Reduction
 */
export class LTSA extends DR {
    /**
     * Local Tangent Space Alignment
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLTSA>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                neighbors: -Infinity,
                d: 2,
                metric: euclidean,
                seed: 1212,
                knn: null,
                eig_args: {},
            },
            parameters,
        );
        if (this.parameter("neighbors") === -Infinity) {
            this.parameter("neighbors", Math.min(Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        }
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
        }

        const d = /** @type {number} */ (this.parameter("d"));
        if (this._D <= d) {
            throw new Error(
                `Dimensionality of X (D = ${this._D}) must be greater than the required dimensionality of the result (d = ${d})!`,
            );
        }
    }

    /**
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {Generator<T, T, void>} A generator yielding the intermediate steps of the projection.
     */
    *generator() {
        yield this.transform();
        return this.projection;
    }

    /**
     * Transforms the inputdata `X` to dimenionality `d`.
     *
     * @returns {T}
     */
    transform() {
        const X = this.X;
        const [rows, D] = X.shape;
        const neighbors = /** @type {number} */ (this.parameter("neighbors"));
        const d = /** @type {number} */ (this.parameter("d"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const metric = /** @type {typeof euclidean} */ (this.parameter("metric"));
        // 1.1 determine k nearest neighbors
        const knn =
            /** @type {KNN<number[] | Float64Array, any> | null} */ (this.parameter("knn")) ??
            spatial_tree(X.to2dArray(), { metric, seed: /** @type {number} */ (this.parameter("seed")) });
        const nN = [];
        for (let row = 0; row < rows; ++row) {
            const found = knn.search_by_index(row, neighbors);
            if (found.length < neighbors) {
                throw new Error(
                    `The KNN index returned only ${found.length} of the ${neighbors} requested neighbors for element ${row}. ` +
                        `Lower "neighbors", or use an index that recalls more candidates.`,
                );
            }
            nN.push(found);
        }
        // center matrix
        const O = new Matrix(D, D, "center");
        const B = new Matrix(rows, rows, 0);

        for (let row = 0; row < rows; ++row) {
            // 1.2 compute the d largest eigenvectors of the correlation matrix
            const I_i = [row, ...nN[row].map((n) => n.index)];
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

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLTSA>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new LTSA(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLTSA>} parameters
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new LTSA(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLTSA>} parameters
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new LTSA(X, parameters);
        return dr.transform_async();
    }
}
