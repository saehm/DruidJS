/**
 * Computes the Cosine distance between vector {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias cosine
 * @param {Array} a 
 * @param {Array} b 
 * @returns {float64} The Cosine distance between vector {@link a} and {@link b}.  
 */
export default function(a, b) {
    if (a.length !== b.length) return undefined;
    let n = a.length;
    let sum = 0;
    let sum_a = 0;
    let sum_b = 0;
    for (let i = 0; i < n; ++i) {
        sum += (a[i] * b[i])
        sum_a += (a[i] * a[i])
        sum_b += (b[i] * b[i])
    }
    return sum / ((Math.sqrt(sum_a) * Math.sqrt(sum_b)));
}