/**
 * Computes the canberra distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias canberra
 * @param {Number[]} a 
 * @param {Number[]} b 
 * @returns {Number} the canberra distance between <code>a</code> and <code>b</code>.
 * @see {@link https://en.wikipedia.org/wiki/Canberra_distance}
 */
export default function(a, b) {
    if (a.length !== b.length) return undefined;
    const n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += (Math.abs(a[i] - b[i]) / (Math.abs(a[i]) + Math.abs(b[i])))
    }
    return sum;
}