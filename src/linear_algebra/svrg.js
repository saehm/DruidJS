import { Matrix, norm } from "../matrix/index";
import { Randomizer } from "../util/randomizer";

/**
 * Computes the eigenvector of {@link X} with an accelerated stochastic power iteration algorithm.
 * @memberof module:linear_algebra 
 * @alias svrg
 * @see {@link https://arxiv.org/abs/1707.02670}
 * @param {Matrix} data - the data matrix
 * @param {Matrix} x - Initial Point as 1 times cols Matrix
 * @param {number} beta - momentum parameter
 * @param {number} epoch - number of epochs
 * @param {number} m - epoch length
 * @param {number} s - mini-batch size
 * @param {number} seed - seed for the random number generator
 */
export default function(data, x, beta, epoch=20, m=10, s=1, seed) {
    let [n, d] = data.shape;
    const randomizer = new Randomizer(seed);
    x = new Matrix(d, 1, () => randomizer.random)
    x = x.divide(norm(x));
    let x0 = x.clone();
    let A = data.transDot(data).divide(n);
    let x_tilde = x.clone();
    
    for (let t = 0; t < epoch; ++t) {
        const gx = A.dot(x_tilde);
        for (let i = 0; i < m; ++i) {
            const ang = x.transDot(x_tilde).entry(0, 0);
            const sample = Matrix.from(Randomizer.choice(data, s));
            const sampleT_dot_sample = sample.transDot(sample)
            const x_tmp = x.clone();
            const XTXx = sampleT_dot_sample
                    .dot(x.divide(s));
            const XTXx_tilde = sampleT_dot_sample
                    .dot(x_tilde.mult(ang / s));
            x = XTXx.sub(XTXx_tilde)
                    .add(gx.mult(ang).sub(x0.mult(beta)));
            x0 = x_tmp;
            const x_norm = norm(x);
            x = x.divide(x_norm);
            x0 = x0.divide(x_norm);        
        }  
        x_tilde = x.clone();
    }
    return x;

}