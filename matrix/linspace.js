/**
 * @memberof module:matrix
 * Creates an Array containing {@link number} numbers from {@link start} to {@link end}. If <code>{@link number} = null</null>
 * @param {Number} start
 * @param {Number} end
 * @param {Number} [number = null]
 * @returns {Array}
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
