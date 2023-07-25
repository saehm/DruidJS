import { distance_matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";

/**
 * Computes the k-nearest neighbors of each row of {@link A}.
 * @memberof module:matrix
 * @alias k_nearest_neigbhors
 * @param {Matrix} A - Either the data matrix, or a distance matrix.
 * @param {Number} k - The number of neighbors to compute.
 * @param {Function|"precomputed"} [metric=euclidean]
 * @returns {Array<Object>} -
 */
export default function (A, k, metric = euclidean) {
    const rows = A.shape[0];
    let D = metric == "precomputed" ? A : distance_matrix(A, metric);
    let nN = new Array(rows);
    for (let row = 0; row < rows; ++row) {
        nN[row] = Array.from(D.row(row))
            .map((distance, col) => {
                return {
                    i: row,
                    j: col,
                    distance: distance,
                };
            })
            .sort((a, b) => a.distance - b.distance)
            .slice(1, k + 1);
    }
    return nN;
}
