/**
 * Numerical stable summation with the Neumair summation algorithm.
 *
 * Deliberately not WASM accelerated: the compensation term makes each step depend on the previous
 * one, so the kernel cannot vectorise, and it measured slower than this loop at every input size —
 * the argument copy is pure overhead. A `neumair_sum_f64` kernel existed on that basis alone, never
 * called by anything; it and its benchmark row were removed once the measurement had been made.
 *
 * @category Numerical
 * @param {number[] | Float64Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements}
 */
export function neumair_sum(summands) {
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
    // On overflow the compensation is `Infinity - Infinity`, so adding it turns a correct `Infinity`
    // into `NaN`. There is nothing to compensate once the running sum has left the representable
    // range, so it is returned as is.
    return Number.isFinite(sum) ? sum + compensation : sum;
}
