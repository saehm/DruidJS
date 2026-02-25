/**
 * Returns maximum in Array `values`.
 *
 * @category Utils
 * @param {Iterable<number | null>} values
 * @returns {number}
 */
export function min(values) {
    let min = Infinity;
    for (const value of values) {
        if (value !== null && min > value) {
            min = value;
        }
    }
    return min;
}
