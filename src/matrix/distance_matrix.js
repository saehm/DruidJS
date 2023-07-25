import { euclidean } from "../metrics/index.js";
import { Matrix } from "./index.js";

/**
 * Computes the distance matrix of datamatrix {@link A}.
 * @memberof module:matrix
 * @alias distance_matrix
 * @param {Matrix} A - Matrix.
 * @param {Function} [metric=euclidean] - The diistance metric.
 * @returns {Matrix} D - The distance matrix of {@link A}.
 */
export default function (A, metric = euclidean) {
    let n = A.shape[0];
    const D = new Matrix(n, n);
    for (let i = 0; i < n; ++i) {
        const A_i = A.row(i);
        for (let j = i + 1; j < n; ++j) {
            const dist = metric(A_i, A.row(j));
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
        }
    }
    return D;
}
