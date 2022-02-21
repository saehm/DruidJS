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
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @returns {PCA}
     */
    constructor(X, d = 2) {
        super(X, d);
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimensionality {@link d}. If parameter {@link A} is given, then project {@link A} with the principal components of {@link X}.
     * @param {null|Matrix|Array} [A = null] - If given, the data to project.
     * @returns {Matrix|Array} - The projected data.
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
        const X = this.X;
        const means = Matrix.from(X.meanCols);
        const X_cent = X.sub(means);
        const C = X_cent.transpose().dot(X_cent);
        const { eigenvectors: V } = simultaneous_poweriteration(C, this._d);
        this.V = Matrix.from(V).transpose();
        return this.V;
    }
}
