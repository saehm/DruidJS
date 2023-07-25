/**
 * Computes the Sokal-Michener distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias sokal_michener
 * @param {Number[]} a 
 * @param {Number[]} b 
 * @returns {Number} the Sokal-Michener distance between <code>a</code> and <code>b</code>.  
 */
export default function(a, b) {
    if (a.length != b.length) return undefined
    const n = a.length;
    let num_not_equal = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] != 0;
        const y = b[i] != 0;
        num_not_equal += x != y;
    }
    return (2 * num_not_equal) / (n + num_not_equal);
}