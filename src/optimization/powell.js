/**
 * @template {Float64Array | number[]} T
 * @category Optimization
 * @param {(d: T) => number} f
 * @param {T} x0
 * @param {number} [max_iter=300] Default is `300`
 * @returns {T}
 * @see http://optimization-js.github.io/optimization-js/optimization.js.html#line438
 */
export function powell(f, x0, max_iter = 300) {
    const epsilon = 1e-2;
    const n = x0.length;
    let alpha = 1e-3;
    let pfx = 10000;
    const x = /** @type {T} */ (x0.slice());
    let fx = f(x);
    let convergence = false;

    while (max_iter-- >= 0 && !convergence) {
        convergence = true;
        for (let i = 0; i < n; ++i) {
            x[i] += 1e-6;
            const fxi = f(x);
            x[i] -= 1e-6;
            const dx = (fxi - fx) / 1e-6;
            if (Math.abs(dx) > epsilon) {
                convergence = false;
            }
            x[i] -= alpha * dx;
            fx = f(x);
        }
        alpha *= pfx >= fx ? 1.05 : 0.4;
        pfx = fx;
    }
    return x;
}
