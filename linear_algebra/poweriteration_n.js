import { norm } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";
import { neumair_sum } from "../numerical/index.js";
/**
 * @typedef {Eigenpair} Eigenpair
 * @property {Array} Eigenvalues - Array of Eigenvalues
 * @property {Array[]} Eigenvectors - Array of Eigenvectors 
 */


/**
 * Computes the {@link n} biggest Eigenpair of the Matrix {@link data}.
 * @memberof module:linear_algebra
 * @alias poweriteration_n
 * @param {Matrix} data - the data matrix
 * @param {int} n - Number of Eigenvalues / Eigenvectors
 * @param {Matrix} x - Initial Point as 1 times cols Matrix
 * @param {number} beta - momentum parameter
 * @param {number} max_iter - maximum number of iterations
 * @param {number} seed - seed for the random number generator
 * @returns {Eigenpair} The {@link n} Eigenpairs.
 */
export default function(data, n, x0, beta, max_iter=100, seed) {
    const randomizer = new Randomizer(seed);
    const N = data.shape[0];
    //let b = new Matrix(N, n, () => randomizer.random);
    let b = [];

    if (x0 == null) {
        x0 = new Array(n)//new Matrix(N, n, () => randomizer.random)
        
        for (let i = 0; i < n; ++i) {
            x0[i] = new Float64Array(N);
            b[i] = new Float64Array(N);
            for (let j = 0; j < N; ++j) {
                const value = randomizer.random
                x0[i][j] = value;
                b[i][j] = value;
            }
            let x0_i_norm = norm(x0[i]);
            x0[i] = x0[i].map(x => x / x0_i_norm);
        }
        //x0 = Matrix.from(x0).T;
        //b = Matrix.from(b).T;
    }
    //x0 = x0.divide(norm(x0));
    for (let k = 0; k < n; ++k) {
        let bk = b[k];
        for (let s = 0; s < max_iter; ++s) {
            // Orthogonalize vector
            for (let l = 0; l < k; ++l) {
                const row = b[l]
                const d = neumair_sum((new Float64Array(N)).map((_, i) => bk[i] * row[i]));
                for (let i = 0; i < N; ++i) {
                    bk[i] = bk[i] - (d * row[i]);
                }
            }
            let tmp = data.dot(bk);
            const tmp_norm = norm(tmp)
            x0[k] = tmp.map(t => t / tmp_norm);
            if (neumair_sum((new Float64Array(N)).map((_, i) => tmp[i] * bk[i])) > (1 - 1e-12)) {
                break;
            }
            [bk, tmp] = [tmp, bk]
        }
    }

    return {
        "eigenvalues": b,
        "eigenvectors": x0
    }
}