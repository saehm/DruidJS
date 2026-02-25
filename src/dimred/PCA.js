import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";
import { DR } from "./DR.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersPCA} from "./index.js" */
/** @import {EigenArgs} from "../linear_algebra/index.js" */

/**
 * Principal Component Analysis (PCA)
 *
 * A linear dimensionality reduction technique that identifies the axes (principal components)
 * along which the variance of the data is maximized.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersPCA>
 * @category Dimensionality Reduction
 * @see {@link MDS} for another linear alternative
 *
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 *
 * const X = [[1, 2], [3, 4], [5, 6]];
 * const pca = new druid.PCA(X, { d: 2 });
 * const Y = pca.transform();
 * // [[x1, y1], [x2, y2], [x3, y3]]
 */
export class PCA extends DR {
    /**
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersPCA>} [parameters] - Object containing parameterization of the DR method.
     */
    constructor(X, parameters = {}) {
        super(X, { d: 2, seed: 1212, eig_args: {} }, parameters);
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
     * @returns {T} - The projected data.
     */
    transform() {
        const V = this.principal_components();
        const X = this.X;
        this.Y = X.dot(V);
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} parameters
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new PCA(X, parameters);
        return dr.transform();
    }

    /**
     * Computes the `d` principal components of Matrix `X`.
     *
     * @returns {Matrix}
     */
    principal_components() {
        if (this.V) {
            return this.V;
        }
        const d = /** @type {number} */ (this.parameter("d"));
        const eig_args = /** @type {Partial<EigenArgs>} */ (this.parameter("eig_args"));
        const X = this.X;
        const X_cent = X.sub(X.meanCols());
        const C = X_cent.transDot(X_cent);
        const { eigenvectors: V } = simultaneous_poweriteration(C, d, eig_args);
        this.V = Matrix.from(V).transpose();
        return this.V;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} parameters
     * @returns {Matrix}
     */
    static principal_components(X, parameters) {
        const dr = new PCA(X, parameters);
        return dr.principal_components();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new PCA(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersPCA>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new PCA(X, parameters);
        return dr.transform_async();
    }
}
