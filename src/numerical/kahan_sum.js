/**
 * Numerical stable summation with the Kahan summation algorithm.
 * @memberof module:numerical
 * @alias kahan_sum
 * @param {Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm}
 */
export default function (summands) {
    let n = summands.length;
    let sum = 0;
    let compensation = 0;
    let y, t;

    for (let i = 0; i < n; ++i) {
        y = summands[i] - compensation;
        t = sum + y;
        compensation = t - sum - y;
        sum = t;
    }
    return sum;
}
