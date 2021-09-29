import { simultaneous_poweriteration} from "../linear_algebra/index";
import { Matrix } from "../matrix/index";
import { DR } from "./DR.js";

/**
 * @class
 * @alias PCA
 * @augments DR
 */
export class PCA extends DR{
    /**
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias PCA 
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @returns {PCA}
     */
    constructor(X, d=2) {
        super(X, d);
        return this;
    }

    /**
     * Returns the principal components of the projection
     */
    get principalCmp() {
        return this.V;
    }

    /**
     * projects the data if the principal components have already been calculated
     * @param {Matrix|Array<Array<Number>>} data - the data to be projected
     */
    project(data) {
        if (Array.isArray(data)) {
            data = Matrix.from(data);
        else if (!(data instanceof Matrix)) {
            throw "no valid type for X";
        }
        return data.dot(this.principalCmp);
    }
    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        let X = this.X;
        let D = X.shape[1];
        let O = new Matrix(D, D, "center");
        let X_cent = X.dot(O);

        let C = X_cent.transpose().dot(X_cent)
        let { eigenvectors: V } = simultaneous_poweriteration(C, this._d)
        this.V = Matrix.from(V).transpose();
        this.Y = X.dot(this.V)
        return this.projection;
    }
} 
