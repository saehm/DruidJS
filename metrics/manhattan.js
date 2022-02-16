/**
 * Computes the manhattan distance (l<sub>1</sub>) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias manhattan
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the manhattan distance between {@link a} and {@link b}.
 */ 
export default function (a, b) {
    if (a.length != b.length) return undefined;
    let n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += Math.abs(a[i] - b[i]);
    }
    return sum;
}
