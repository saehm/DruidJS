/**
 * Computes the yule distance between `a` and `b`.
 *
 * @category Metrics
 * @param {number[] | Float64Array} a
 * @param {number[] | Float64Array} b
 * @returns {number} The yule distance between `a` and `b`.
 */
export function yule(a, b) {
    if (a.length !== b.length) throw new Error("Vector a and b needs to be of the same length!");
    const n = a.length;
    let num_true_true = 0;
    let num_true_false = 0;
    let num_false_true = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] !== 0;
        const y = b[i] !== 0;
        num_true_true += x && y ? 1 : 0;
        num_true_false += x && !y ? 1 : 0;
        num_false_true += !x && y ? 1 : 0;
    }
    const num_false_false = n - num_true_true - num_true_false - num_false_true;
    return num_true_false === 0 || num_false_true === 0
        ? 0
        : (2 * num_true_false * num_false_true) / (num_true_true * num_false_false + num_true_false * num_false_true);
}
