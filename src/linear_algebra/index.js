/** @module linear_algebra */

/** @import {Randomizer} from "../util/index.js" */
/**
 * @callback QRDecomposition
 * @param {import("../matrix/index.js").Matrix} A
 * @returns {{ R: import("../matrix/index.js").Matrix, Q: import("../matrix/index.js").Matrix }}
 */

/**
 * @typedef {Object} EigenArgs
 * @property {number} [max_iterations=100] - The number of maxiumum iterations the algorithm should run. Default is `100`
 * @property {number | Randomizer} [seed=1212] - The seed value or a randomizer used in the algorithm. Default is `1212`
 * @property {QRDecomposition} [qr] - The QR technique to use. Default is `qr_gramschmidt`
 * @property {number} [tol=1e-8] - Tolerated error for stopping criteria. Default is `1e-8`
 */

export { inner_product } from "./inner_product.js";
export { qr } from "./qr.js";
export { qr_householder } from "./qr_householder.js";
//export { default as qr_givens } from "./qr_givens.js";
export { simultaneous_poweriteration } from "./simultaneous_poweriteration.js";
