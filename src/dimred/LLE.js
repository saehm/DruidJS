import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { k_nearest_neighbors, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { neumair_sum } from "../numerical/index.js";
import { DR } from "./DR.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersLLE} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */

/**
 * Locally Linear Embedding (LLE)
 *
 * A nonlinear dimensionality reduction technique that preserves local
 * linear relationships between points. It represents each point as a linear
 * combination of its neighbors.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersLLE>
 * @category Dimensionality Reduction
 * @see {@link ISOMAP} for another nonlinear alternative
 */
export class LLE extends DR {
    /**
     * Locally Linear Embedding.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLLE>} parameters - Object containing parameterization of the DR method.
     * @see {@link https://doi.org/10.1126/science.290.5500.2323}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                neighbors: -Infinity,
                d: 2,
                metric: euclidean,
                seed: 1212,
                eig_args: {},
            },
            parameters,
        );
        if (this._parameters.neighbors === -Infinity) {
            this.parameter("neighbors", Math.min(Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        }

        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        if (!Object.hasOwn(eig_args, "seed")) {
            eig_args.seed = this._randomizer;
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
     * Transforms the inputdata `X` to dimensionality `d`.
     *
     * @returns {T}
     */
    transform() {
        const X = this.X;
        const rows = this._N;
        const cols = this._D;
        const neighbors = /** @type {number} */ (this.parameter("neighbors"));
        const d = /** @type {number} */ (this.parameter("d"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const metric = /** @type {typeof euclidean} */ (this.parameter("metric"));
        const nN = k_nearest_neighbors(X, neighbors, metric);
        const O = new Matrix(neighbors, 1, 1);
        const W = new Matrix(rows, rows);

        for (let row = 0; row < rows; ++row) {
            const nN_row = nN[row];
            const Z = new Matrix(neighbors, cols, (i, j) => X.entry(nN_row[i].j, j) - X.entry(row, j));
            const C = Z.dotTrans(Z);
            if (neighbors > cols) {
                const C_trace = neumair_sum(C.diag()) / 1000;
                for (let j = 0; j < neighbors; ++j) {
                    C.add_entry(j, j, C_trace);
                }
            }
            // reconstruct;
            let w = Matrix.solve_CG(C, O, this._randomizer);
            w = w.divide(w.sum());
            for (let j = 0; j < neighbors; ++j) {
                W.set_entry(row, nN_row[j].j, w.entry(j, 0));
            }
        }
        // comp embedding
        const I = new Matrix(rows, rows, "identity");
        const IW = I.sub(W);
        const M = IW.transDot(IW);

        // M is symmetric positive semi-definite. Smallest eigenvalue is 0 (ones vector).
        // To find smallest eigenvalues of M, we can find largest of (C*I - M)
        // Upper bound for max eigenvalue: Frobenius norm or sum of absolute values
        const C = M.mean() * rows * 2; // Safe upper bound for a sparse-ish M in LLE
        const CI_M = new Matrix(rows, rows, (i, j) => (i === j ? C : 0) - M.entry(i, j));

        const { eigenvectors: V } = simultaneous_poweriteration(CI_M, d + 1, eig_args);
        // Skip the first eigenvector (the ones vector corresponding to eigenvalue C)
        this.Y = Matrix.from(V.slice(1, 1 + d)).T;

        // return embedding
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLLE>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new LLE(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLLE>} parameters
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new LLE(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLLE>} parameters
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new LLE(X, parameters);
        return dr.transform_async();
    }
}
