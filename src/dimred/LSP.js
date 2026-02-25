import { KMedoids } from "../clustering/index.js";
import { BallTree } from "../knn/index.js";
import { Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { DR } from "./DR.js";
import { MDS } from "./MDS.js";

/** @import {InputType} from "../index.js" */
/** @import {ParametersLSP} from "./index.js" */

/**
 * Least Square Projection (LSP)
 *
 * A dimensionality reduction technique that uses a small set of control points
 * (projected with MDS) to define the projection for the rest of the data
 * using a Laplacian-based optimization.
 *
 * @class
 * @template {InputType} T
 * @extends DR<T, ParametersLSP>
 * @category Dimensionality Reduction
 */
export class LSP extends DR {
    /**
     * Least Squares Projection.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Partial<ParametersLSP>} [parameters] - Object containing parameterization of the DR method.
     * @see {@link https://ieeexplore.ieee.org/document/4378370}
     */
    constructor(X, parameters) {
        super(
            X,
            {
                neighbors: -Infinity,
                control_points: -Infinity,
                d: 2,
                metric: euclidean,
                seed: 1212,
            },
            parameters,
        );
        if (this.parameter("neighbors") === -Infinity) {
            this.parameter("neighbors", Math.min(Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        }
        if (this.parameter("control_points") === -Infinity) {
            this.parameter("control_points", Math.min(Math.ceil(Math.sqrt(this._N)), this._N - 1));
        }
        this._is_initialized = false;
    }

    /**
     * @returns {LSP<T>}
     */
    //	init(DR = MDS, DR_parameters = {}, KNN = BallTree) {
    init() {
        const DR = MDS;
        let DR_parameters = {};
        const KNN = BallTree;
        if (this._is_initialized) return this;
        const X = this.X;
        const N = this._N;
        const K = /** @type {number} */ (this.parameter("neighbors"));
        const d = /** @type {number} */ (this.parameter("d"));
        const seed = /** @type {number} */ (this.parameter("seed"));
        const metric = /** @type {typeof euclidean} */ (this.parameter("metric"));
        DR_parameters = Object.assign({ d, metric, seed }, DR_parameters);
        const nc = /** @type {number} */ (this.parameter("control_points"));
        const control_points = new KMedoids(X, { K: nc, metric }).get_medoids();
        const C = new Matrix(nc, N, "zeros");
        control_points.forEach((c_i, i) => {
            C.set_entry(i, c_i, 1);
        });

        const control_points_matrix = Matrix.from(control_points.map((c_i) => X.row(c_i)));
        const Y_C = new DR(control_points_matrix, DR_parameters).transform();

        const XA = X.to2dArray();
        const knn = new KNN(XA, { metric, seed });
        const L = new Matrix(N, N, "I");
        const alpha = -1 / K;
        XA.forEach((x_i, i) => {
            for (const { index: j } of knn.search(x_i, K)) {
                if (i === j) continue;
                L.set_entry(i, j, alpha);
            }
        });
        const A = L.concat(C, "vertical");

        const z = new Matrix(N, d, "zeros");
        const b = z.concat(Y_C, "vertical");

        this._A = A;
        this._b = b;
        this._is_initialized = true;
        return this;
    }

    /**
     * Computes the projection.
     *
     * @returns {T} Returns the projection.
     */
    transform() {
        this.check_init();
        const A = this._A;
        const b = this._b;

        if (!A || !b) throw new Error("Call init() first!");
        const ATA = A.transDot(A);
        const ATb = A.transDot(b);
        this.Y = Matrix.solve_CG(ATA, ATb, this._randomizer);
        return this.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLSP>} [parameters]
     * @returns {T}
     */
    static transform(X, parameters) {
        const dr = new LSP(X, parameters);
        return dr.transform();
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLSP>} [parameters]
     * @returns {Generator<T, T, void>}
     */
    static *generator(X, parameters) {
        const dr = new LSP(X, parameters);
        yield* dr.generator();
        return dr.projection;
    }

    /**
     * @template {InputType} T
     * @param {T} X
     * @param {Partial<ParametersLSP>} [parameters]
     * @returns {Promise<T>}
     */
    static async transform_async(X, parameters) {
        const dr = new LSP(X, parameters);
        return dr.transform_async();
    }
}
