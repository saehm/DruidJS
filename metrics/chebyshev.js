/**
 * Computes the Chebyshev distance between vector {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias chebyshev
 * @param {Array} a 
 * @param {Array} b 
 * @returns {float64} the Chebyshev distance between vector {@link a} and {@link b}.  
 */
export default function(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length
    let res = [];
    for (let i = 0; i < n; ++i) {
        res.push(Math.abs(a[i] - b[i]))
    }
    return Math.max(...res)
}