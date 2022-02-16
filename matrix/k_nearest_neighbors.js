import { distance_matrix as dmatrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";

/**
 *
 * @param {*} A
 * @param {*} k
 * @param {*} distance_matrix
 * @param {*} metric
 */
export default function (A, k, distance_matrix = null, metric = euclidean) {
    const rows = A.shape[0];
    let D = distance_matrix ?? dmatrix(A, metric);
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
