import { simultaneous_poweriteration } from "../linear_algebra/index.js";
import { Matrix } from "../matrix/index.js";

/** @import { EigenArgs } from "../linear_algebra/index.js" */

/**
 * Classical (Torgerson) MDS of a distance matrix: `Y = V Λ^½` for the `d` largest eigenvalues `λ`
 * of the double-centered matrix `B = -½ J D² J`.
 *
 * Distances that are not Euclidean (precomputed dissimilarities, Isomap's geodesic distances) give
 * `B` negative eigenvalues. Power iteration finds the eigenvalues largest in magnitude, so when a
 * negative one turns up, it is the most negative of all: shifting `B` by it makes `B + cI` positive
 * semidefinite, and a second run finds the largest eigenvalues by value. Axes with `λ ≤ 0` get no
 * extent. Euclidean distances never need the second run.
 *
 * @private
 * @param {Matrix} D - The `n × n` distance matrix.
 * @param {number} d - The target dimensionality.
 * @param {Partial<EigenArgs>} [eig_args] - Parameters for {@link simultaneous_poweriteration}.
 * @returns {Matrix} The `n × d` embedding.
 */
export function classical_mds(D, d, eig_args = {}) {
    const n = D.shape[0];
    const D_sq = new Matrix(n, n, (i, j) => D.entry(i, j) ** 2);
    const ai_ = D_sq.meanCols();
    const a_j = D_sq.meanRows();
    const a__ = D_sq.mean();
    const B = new Matrix(n, n, (i, j) => -0.5 * (D_sq.entry(i, j) - ai_[i] - a_j[j] + a__));

    let { eigenvalues, eigenvectors } = simultaneous_poweriteration(B, d, eig_args);
    const most_negative = Math.min(...eigenvalues);
    if (most_negative < 0) {
        const shift = -most_negative;
        const B_shifted = new Matrix(n, n, (i, j) => B.entry(i, j) + (i === j ? shift : 0));
        ({ eigenvalues, eigenvectors } = simultaneous_poweriteration(B_shifted, d, eig_args));
        eigenvalues = eigenvalues.map((lambda) => lambda - shift);
    }

    const scales = Array.from(eigenvalues, (lambda) => Math.sqrt(Math.max(lambda, 0)));
    return new Matrix(n, d, (i, k) => eigenvectors[k][i] * scales[k]);
}
