/**
 * Computes the yule distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias yule
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the yule distance between <code>a</code> and <code>b</code>.
 */
export default function (a, b) {
    if (a.length != b.length) return undefined;
    const n = a.length;
    let num_true_true = 0;
    let num_true_false = 0;
    let num_false_true = 0;
    for (let i = 0; i < n; ++i) {
        const x = a[i] != 0;
        const y = b[i] != 0;
        num_true_true += x && y;
        num_true_false += x && !y;
        num_false_true += !x && x;
    }
    const num_false_false = n - num_true_true - num_true_false - num_false_true;
    return num_true_false == 0 || num_false_true == 0 ? 0 : (2 * num_true_false * num_false_true) / (num_true_true * num_false_false + num_true_false * num_false_true);
}
