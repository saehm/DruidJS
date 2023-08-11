import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";
import { DR } from "./DR.js";

/**
 * @class
 * @alias PCA
 * @augments DR
 */
export class PCA extends DR {
    /**
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias PCA
     * @param {Matrix|number[][]} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {object} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @returns {PCA}
     */
    constructor(X, parameters) {
        super(X, { d: 2, seed: 1212, eig_args: {} }, parameters);
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimensionality {@link d}. If parameter {@link A} is given, then project {@link A} with the principal components of {@link X}.
     * @param {null|Matrix|number[][]} [A = null] - If given, the data to project.
     * @returns {Matrix|number[][]} - The projected data.
     */
    transform(A = null) {
        const V = this.principal_components();
        if (A == null) {
            const X = this.X;
            this.Y = X.dot(V);
            return this.projection;
        } else if (Array.isArray(A)) {
            return Matrix.from(A).dot(V).asArray;
        } else if (A instanceof Matrix) {
            return A.dot(V);
        } else {
            throw new Error("No valid type for A!");
        }
    }

    /**
     * Computes the {@link d} principal components of Matrix {@link X}.
     * @returns {Matrix}
     */
    principal_components() {
        if (this.V) {
            return this.V;
        }
        const { d, eig_args } = this._parameters;
        const X = this.X;
        const X_cent = X.sub(X.meanCols);
        const C = X_cent.transDot(X_cent);
        const { eigenvectors: V } = simultaneous_poweriteration(C, d, eig_args);
        this.V = Matrix.from(V).transpose();
        return this.V;
    }

    static principal_components(X, parameters) {
        const dr = new this(X, parameters);
        return dr.principal_components();
    }
}
