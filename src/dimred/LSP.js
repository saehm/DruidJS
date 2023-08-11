import { Matrix } from "../matrix/index.js";
import { DR } from "./DR.js";
import { MDS } from "./MDS.js";
import { KMedoids } from "../clustering/index.js";
import { euclidean } from "../metrics/index.js";
import { BallTree } from "../knn/index.js";
/**
 * @class
 * @alias LSP
 * @extends DR
 */
export class LSP extends DR {
    /**
     * Least Squares Projection.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LSP
     * @param {Matrix} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} [parameters.neighbors = Math.max(Math.floor(N / 10), 2)] - number of neighbors to consider.
     * @param {number} [parameters.control_points = Math.ceil(Math.sqrt(N))] - number of controlpoints
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {LSP}
     * @see {@link https://ieeexplore.ieee.org/document/4378370}
     * @todo accept precomputed distance matrix.
     */
    constructor(X, parameters) {
        super(X, { neighbors: undefined, control_points: undefined, d: 2, metric: euclidean, seed: 1212 }, parameters);
        this.parameter("neighbors", Math.min(this._parameters.neighbors ?? Math.max(Math.floor(this._N / 10), 2), this._N - 1));
        this.parameter("control_points", Math.min(this._parameters.control_points ?? Math.ceil(Math.sqrt(this._N)), this._N - 1));
        this._is_initialized = false;
        return this;
    }

    /**
     *
     * @param {DR} DR - method used for position control points.
     * @param {object} DR_parameters - Object containing parameters for the DR method which projects the control points
     * @returns {LSP}
     */
    init(DR = MDS, DR_parameters = {}, KNN = BallTree) {
        if (this._is_initialized) return this;
        const X = this.X;
        const N = this._N;
        const K = this.parameter("neighbors");
        const d = this.parameter("d");
        const seed = this.parameter("seed");
        const metric = this.parameter("metric");
        DR_parameters = Object.assign({ d, metric, seed }, DR_parameters);
        const nc = this.parameter("control_points");
        const control_points = new KMedoids(X, nc, null, metric).get_clusters().medoids;
        const C = new Matrix(nc, N, "zeros");
        control_points.forEach((c_i, i) => {
            C.set_entry(i, c_i, 1);
        });
        const Y_C = new DR(Matrix.from(control_points.map((c_i) => X.row(c_i))), DR_parameters).transform();

        const XA = X.to2dArray;
        const knn = new KNN(XA, metric);
        const L = new Matrix(N, N, "I");
        const alpha = -1 / K;
        XA.forEach((x_i, i) => {
            for (const { index: j } of knn.search(x_i, K).iterate()) {
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
     * @returns {Matrix} Returns the projection.
     */
    transform() {
        this.check_init();
        const A = this._A;
        const b = this._b;
        const ATA = A.transDot(A);
        const ATb = A.transDot(b);
        this.Y = Matrix.solve_CG(ATA, ATb, this._randomizer);
        return this.projection;
    }
}
