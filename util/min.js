/**
 * Returns maximum in Array {@link values}.
 * @memberof module:utils
 * @alias min
 * @param {Array} values
 * @returns {Number}
 */
export default function (values) {
    let min;
    for (const value of values) {
        if (value != null && (min > value || (min === undefined && value <= value))) {
            min = value;
        }
    }
    return min;
}