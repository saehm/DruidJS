/**
 * @class
 * @alias DR
 * @borrows DR#parameter as DR#para
 * @borrows DR#parameter as DR#p
 */
export class DR {
    /**
     * @static
     * @param  {...any} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Matrix|Array} - The dimensionality reduced dataset.
     */
    static transform(...args: any[]): Matrix | any[];
    /**
     * @static
     * @param  {...any} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Promise} - A promise yielding the dimensionality reduced dataset.
     */
    static transform_async(...args: any[]): Promise<any>;
    /**
     * @static
     * @param  {...any} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Generator} - A generator yielding the intermediate steps of the dimensionality reduction method.
     */
    static generator(...args: any[]): Generator;
    /**
     * Takes the default parameters and seals them, remembers the type of input {@link X}, and initializes the random number generator.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias DR
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed value for the random number generator.
     * @returns {DR}
     */
    constructor(X: Matrix | Array<Array<number>>, default_parameters: any, parameters: {
        d?: number;
        metric?: Function;
        seed?: number;
    });
    _parameters: any;
    _type: string;
    X: Matrix;
    _randomizer: Randomizer;
    _is_initialized: boolean;
    /**
     * Set and get parameters
     * @param {String} name - name of the parameter.
     * @param {any} [value = null] - value of the parameter to set.
     * @returns {DR|any} - On setting a parameter, this function returns the DR object. If <code>value == null</code> then return actual parameter value.
     * @example
     * const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
     * DR.parameter("d"); // returns 3,
     * DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
     */
    parameter(name: string, value?: any): DR | any;
    para(name: any, value?: any): any;
    p(name: any, value?: any): any;
    /**
     * Computes the projection.
     * @returns {Matrix} - Returns the projection.
     */
    transform(): Matrix;
    /**
     * Computes the projection.
     * @returns {Generator} - A generator yielding the intermediate steps of the dimensionality reduction method.
     */
    generator(): Generator;
    /**
     * If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.
     * @returns {DR}
     */
    check_init(): DR;
    /**
     * @returns {Matrix|Array} Returns the projection.
     */
    get projection(): any[] | Matrix;
    /**
     *
     * @param  {...any} args - Arguments the transform method of the respective DR method takes.
     * @returns {Promise} - A promise yielding the dimensionality reduced dataset.
     */
    transform_async(...args: any[]): Promise<any>;
}
import { Matrix } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";
//# sourceMappingURL=DR.d.ts.map