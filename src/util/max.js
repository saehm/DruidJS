/**
 * Returns maximum in Array `values`.
 *
 * @category Utils
 * @param {Iterable<number | null>} values
 * @returns {number}
 */
export function max(values) {
    let max = -Infinity;
    for (const value of values) {
        if (value !== null && max < value) {
            max = value;
        }
    }
    return max;
}
