/**
 * Creates an Array containing `number` numbers from `start` to `end`. If `number = null`.
 *
 * @category Matrix
 * @param {number} start - Start value.
 * @param {number} end - End value.
 * @param {number} [number] - Number of number between `start` and `end`.
 * @returns {number[]} An array with `number` entries, beginning at `start` ending at `end`.
 */
export function linspace(start, end, number) {
    if (number === undefined || number === null) {
        number = Math.max(Math.round(end - start) + 1, 1);
    }
    if (number < 2) {
        return number === 1 ? [start] : [];
    }
    const result = new Array(number);
    number -= 1;
    for (let i = number; i >= 0; --i) {
        result[i] = (i * end + (number - i) * start) / number;
    }
    return result;
}
