//@ts-check

import { distance_matrix, Matrix } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";

/** @import { Metric } from "../metrics/index.js" */

/**
 * Computes the k-nearest neighbors of each row of `A`.
 *
 * @category Matrix
 * @param {Matrix} A - Either the data matrix, or a distance matrix.
 * @param {number} k - The number of neighbors to compute.
 * @param {Metric | "precomputed"} [metric=euclidean] Default is `euclidean`
 * @returns {{ i: number; j: number; distance: number }[][]} The kNN graph.
 */
export function k_nearest_neighbors(A, k, metric = euclidean) {
    A = A instanceof Matrix ? A : Matrix.from(A);
    const rows = A.shape[0];
    const D = metric === "precomputed" ? A : distance_matrix(A, metric);
    /** @type {{ i: number; j: number; distance: number }[][]} */
    const nN = [];
    for (let row = 0; row < rows; ++row) {
        const res = Array.from(D.row(row))
            .map((distance, col) => {
                return {
                    i: row,
                    j: col,
                    distance: distance,
                };
            })
            .sort((a, b) => a.distance - b.distance)
            .slice(1, k + 1);
        nN.push(res);
    }
    return nN;
}
