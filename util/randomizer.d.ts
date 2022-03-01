/**
 * @class
 * @memberof module:utils
 * @alias Randomizer
 */
export class Randomizer {
    /**
     * @static
     * Returns samples from an input Matrix or Array.
     * @param {Matrix|Number[]|Float64Array} A - The input Matrix or Array.
     * @param {Number} n - The number of samples.
     * @param {Number|Randomizer} seed - The seed for the random number generator.
     * @returns {Number[]} - A random selection form {@link A} of {@link n} samples.
     */
    static choice(A: Matrix | number[] | Float64Array, n: number, seed?: number | Randomizer): number[];
    /**
     * Mersenne Twister random number generator.
     * @constructor
     * @param {Number} [_seed=new Date().getTime()] - The seed for the random number generator. If <code>_seed == null</code> then the actual time gets used as seed.
     * @see https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js
     */
    constructor(_seed?: number);
    _N: number;
    _M: number;
    _MATRIX_A: number;
    _UPPER_MASK: number;
    _LOWER_MASK: number;
    _mt: any[];
    _mti: any;
    set seed(arg: number);
    /**
     * Returns the seed of the random number generator.
     * @returns {Number} - The seed.
     */
    get seed(): number;
    _seed: number;
    /**
     * Returns a float between 0 and 1.
     * @returns {Number} - A random number between [0, 1]
     */
    get random(): number;
    /**
     * Returns an integer between 0 and MAX_INTEGER.
     * @returns {Number} - A random integer.
     */
    get random_int(): number;
    /**
     * Returns samples from an input Matrix or Array.
     * @param {Matrix|Number[]|Float64Array} A - The input Matrix or Array.
     * @param {Number} n - The number of samples.
     * @returns {Number[]} - A random selection form {@link A} of {@link n} samples.
     */
    choice(A: Matrix | number[] | Float64Array, n: number): number[];
}
import { Matrix } from "../matrix/index.js";
//# sourceMappingURL=randomizer.d.ts.map