/**
 * Computes the manhattan distance (<code>l<sub>1</sub></code>) between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias manhattan
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the manhattan distance between <code>a</code> and <code>b</code>.
 */ 
export default function (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += Math.abs(a[i] - b[i]);
    }
    return sum;
}
