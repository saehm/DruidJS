/**
 * Returns maximum in Array {@link values}.
 * @memberof module:utils
 * @alias max
 * @param {Array} values 
 * @returns {Number}
 */
export default function (values) {
    let max;
    for (const value of values) {
        if (value != null && (max < value || (max === undefined && value >= value))) {
            max = value;
        }
    }
    return max;
}