import { Matrix, norm } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";

/**
 * 
 * @param {Matrix} data - the data matrix
 * @param {Matrix} x - Initial Point as 1 times cols Matrix
 * @param {number} beta - momentum parameter
 * @param {number} max_iter - maximum number of iterations
 * @param {number} seed - seed for the random number generator
 */
export default function(data, x0, beta, max_iter=20, seed) {
    let randomizer = new Randomizer(seed);
    let [ n, d ] = data.shape;
    let A = data.transDot(data).divide(n)
    if (x0 === null) x0 = new Matrix(d, 1, () => randomizer.random)
    x0 = x0.divide(norm(x0));
    let x = x0.clone()
    for (let i = 0; i < max_iter; ++i) {
        let x_tmp = x.clone();
        x = A.dot(x).sub(x0.mult(beta));
        x0 = x_tmp;
        let z = norm(x);
        x = x.divide(z);
        x0 = x0.divide(z);
    }
    return x;
}