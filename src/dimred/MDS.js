import { distance_matrix, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { classical_mds } from "./classical_mds.js";
import { DR } from "./DR.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersMDS} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */

/**
 * Classical Multidimensional Scaling (MDS)
 *
 * A linear dimensionality reduction technique that seeks to preserve the
 * pairwise distances between points as much as possible in the lower-dimensional
 * space.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersMDS>
 * @category Dimensionality Reduction
 * @see {@link PCA} for another linear alternative
 */
export class MDS extends DR {
    /**
     * Classical MDS.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersMDS>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters = {}) {
        super(X, { d: 2, metric: euclidean, seed: 1212, eig_args: {} }, parameters);
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
        const d = /** @type {number} */ (this.parameter("d"));
        const metric = /** @type {typeof euclidean | "precomputed"} */ (this.parameter("metric"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const A = metric === "precomputed" ? X : distance_matrix(X, metric);

        this._d_X = A;
        this.Y = classical_mds(A, d, eig_args);

        return this.projection;
    }

    /** @returns {number} - The stress of the projection. */
    stress() {
        const N = this.X.shape[0];
        const Y = this.Y;
        const d_X = this._d_X;
        if (!d_X) throw new Error("First transform!");

        const d_Y = new Matrix(N, N, 0);
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
                top_sum += (d_X.entry(i, j) - d_Y.entry(i, j)) ** 2;
                bottom_sum += d_X.entry(i, j) ** 2;
            }
        }
        return Math.sqrt(top_sum / bottom_sum);
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMDS>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new MDS(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMDS>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new MDS(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersMDS>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new MDS(X, parameters);
        return dr.transform_async();
    }
}
