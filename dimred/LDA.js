import { Matrix } from "../matrix/index";
import { euclidean } from "../metrics/index";
import { simultaneous_poweriteration} from "../linear_algebra/index";
import { DR } from "./DR.js";

/**
 * @class
 * @alias LDA
 */
export class LDA extends DR {
    /**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LDA
     * @param {Matrix} X - the high-dimensional data.
     * @param {Array} labels - the label / class of each data point.
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
    constructor(X, labels, d = 2, metric = euclidean, seed=1212) {
        super(X, d, metric, seed);
        super.parameter_list = ["labels"];
        this.parameter("labels", labels);
        return this;
    }

    /**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */
    transform() {
        let X = this.X;
        let [ rows, cols ] = X.shape;
        let labels = this._labels;
        let unique_labels = {};
        let label_id = 0;
        labels.forEach((l, i) => {
            if (l in unique_labels) {
                unique_labels[l].count++;
                unique_labels[l].rows.push(X.row(i));
            } else {
                unique_labels[l] = {
                    "id": label_id++,
                    "count": 1,
                    "rows": [X.row(i)]
                };
            }
        })
        
        // create X_mean and vector means;
        let X_mean = X.mean;
        let V_mean = new Matrix(label_id, cols)
        for (let label in unique_labels) {
            let V = Matrix.from(unique_labels[label].rows);
            let v_mean = V.meanCols;
            for (let j = 0; j < cols; ++j) {
                V_mean.set_entry(unique_labels[label].id, j, v_mean[j]);
            }           
        }
        // scatter_between
        let S_b = new Matrix(cols, cols);
        for (let label in unique_labels) {
            let v = V_mean.row(unique_labels[label].id);
            let m = new Matrix(cols, 1, (j) => v[j] - X_mean);
            let N = unique_labels[label].count;
            S_b = S_b.add(m.dot(m.transpose()).mult(N));
        }

        // scatter_within
        let S_w = new Matrix(cols, cols);
        for (let label in unique_labels) {
            let v = V_mean.row(unique_labels[label].id);
            let m = new Matrix(cols, 1, (j) => v[j])
            let R = unique_labels[label].rows;
            for (let i = 0, n = unique_labels[label].count; i < n; ++i) {
                let row_v = new Matrix(cols, 1, (j,_) => R[i][j] - m.entry(j, 0));
                S_w = S_w.add(row_v.dot(row_v.transpose()));
            }
        }

        let { eigenvectors: V } = simultaneous_poweriteration(S_w.inverse().dot(S_b), this._d)
        V = Matrix.from(V).transpose()
        this.Y = X.dot(V)

        // return embedding
        return this.projection;
    }
}