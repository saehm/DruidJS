/**
 * Numerical stable summation with the Neumair summation algorithm.
 * @memberof module:numerical
 * @alias neumair_sum
 * @param {Number[]} summands - Array of values to sum up.
 * @returns {Number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements}
 */
export default function (summands) {
    const n = summands.length;
    let sum = 0;
    let compensation = 0;

    for (let i = 0; i < n; ++i) {
        const summand = summands[i];
        const t = sum + summand;
        if (Math.abs(sum) >= Math.abs(summand)) {
            compensation += sum - t + summand;
        } else {
            compensation += summand - t + sum;
        }
        sum = t;
    }
    return sum + compensation;
}
