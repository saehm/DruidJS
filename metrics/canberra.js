/**
 * Computes the canberra distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias canberra
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} The canberra distance between {@link a} and {@link b}.
 * @see {@link https://en.wikipedia.org/wiki/Canberra_distance}
 */
export default function(a, b) {
    if (a.length !== b.length) return undefined;
    let n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += (Math.abs(a[i] - b[i]) / (Math.abs(a[i]) + Math.abs(b[i])))
    }
    return sum;
}