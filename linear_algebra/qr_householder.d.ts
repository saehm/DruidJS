/**
 * Computes the QR Decomposition of the Matrix {@link A} with householder transformations.
 * @memberof module:linear_algebra
 * @alias qr_householder
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections}
 * @see {@link http://mlwiki.org/index.php/Householder_Transformation}
 */
export default function _default(A: Matrix): {
    R: Matrix;
    Q: Matrix;
};
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=qr_householder.d.ts.map