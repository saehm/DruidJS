import { euclidean_squared } from "../metrics/index";
/**
 * Computes the euclidean distance (l_2) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias euclidean
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the euclidean distance between {@link a} and {@link b}.  
 */
export default function(a, b) {
    return Math.sqrt(euclidean_squared(a, b));
}