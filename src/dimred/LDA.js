import { Matrix } from "../matrix/index.js";
import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { DR } from "./DR.js";

/**
 * @class
 * @alias LDA
 * @extends DR
 */
export class LDA extends DR {
    /**
     * Linear Discriminant Analysis.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LDA
     * @param {Matrix} X - The high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {any[]} parameters.labels - The labels / classes for each data point.
     * @param {number} [parameters.d = 2] - The dimensionality of the projection.
     * @param {number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {object} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x}
     */
    constructor(X, parameters) {
        super(X, { labels: null, d: 2, seed: 1212, eig_args: {} }, parameters);
        if (!this._parameters.eig_args.hasOwnProperty("seed")) {
            this._parameters.eig_args.seed = this._randomizer;
        }
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        const X = this.X;
        const [rows, cols] = X.shape;
        const { d, labels, eig_args } = this._parameters;
        if (labels === null || labels.length != rows) {
            throw new Error("LDA needs parameter label to every datapoint to work!");
        }
        const unique_labels = {};
        let label_id = 0;
        labels.forEach((l, i) => {
            if (l in unique_labels) {
                unique_labels[l].count++;
                unique_labels[l].rows.push(X.row(i));
            } else {
                unique_labels[l] = {
                    id: label_id++,
                    count: 1,
                    rows: [X.row(i)],
                };
            }
        });

        // create X_mean and vector means;
        const X_mean = X.mean;
        const V_mean = new Matrix(label_id, cols);
        for (const label in unique_labels) {
            const V = Matrix.from(unique_labels[label].rows);
            const v_mean = V.meanCols;
            for (let j = 0; j < cols; ++j) {
                V_mean.set_entry(unique_labels[label].id, j, v_mean[j]);
            }
        }
        // scatter_between
        let S_b = new Matrix(cols, cols);
        for (const label in unique_labels) {
            const v = V_mean.row(unique_labels[label].id);
            const m = new Matrix(cols, 1, (j) => v[j] - X_mean);
            const N = unique_labels[label].count;
            S_b = S_b.add(m.dotTrans(m).mult(N));
        }

        // scatter_within
        let S_w = new Matrix(cols, cols);
        for (const label in unique_labels) {
            const v = V_mean.row(unique_labels[label].id);
            const m = new Matrix(cols, 1, (j) => v[j]);
            const R = unique_labels[label].rows;
            for (let i = 0, n = unique_labels[label].count; i < n; ++i) {
                const row_v = new Matrix(cols, 1, (j, _) => R[i][j] - m.entry(j, 0));
                S_w = S_w.add(row_v.dotTrans(row_v));
            }
        }

        let { eigenvectors: V } = simultaneous_poweriteration(S_w.inverse().dot(S_b), d, eig_args);
        V = Matrix.from(V).transpose();
        this.Y = X.dot(V);

        // return embedding
        return this.projection;
    }
}
