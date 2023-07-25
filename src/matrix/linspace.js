/**
 * Creates an Array containing {@link number} numbers from {@link start} to {@link end}.
 * If <code>{@link number} = null</null>.
 * @memberof module:matrix
 * @alias linspace
 * @param {Number} start - Start value.
 * @param {Number} end - End value.
 * @param {Number} [number = null] - Number of number between {@link start} and {@link end}.
 * @returns {Array} - An array with {@link number} entries, beginning at {@link start} ending at {@link end}.
 */
export default function (start, end, number = null) {
    if (!number) {
        number = Math.max(Math.round(end - start) + 1, 1);
    }
    if (number < 2) {
        return number === 1 ? [start] : [];
    }
    let result = new Array(number);
    number -= 1;
    for (let i = number; i >= 0; --i) {
        result[i] = (i * end + (number - i) * start) / number;
    }
    return result;
}
