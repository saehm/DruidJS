import { Matrix } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";

/**
 * @class
 * @alias DR
 * @borrows DR#parameter as DR#para
 * @borrows DR#parameter as DR#p
 */
export class DR {
    /**
     * Takes the default parameters and seals them, remembers the type of input {@link X}, and initializes the random number generator.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias DR
     * @param {Matrix|number[][]} X - the high-dimensional data.
     * @param {object} parameters - Object containing parameterization of the DR method.
     * @param {number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {number} [parameters.seed = 1212] - the seed value for the random number generator.
     * @returns {DR}
     */
    constructor(X, default_parameters, parameters) {
        this._parameters = Object.assign(Object.seal(default_parameters), parameters);
        if (Array.isArray(X)) {
            this._type = "array";
            this.X = Matrix.from(X);
        } else if (X instanceof Matrix) {
            this._type = "matrix";
            this.X = X;
        } else {
            throw new Error("No valid type for X!");
        }
        [this._N, this._D] = this.X.shape;
        this._randomizer = new Randomizer(this._parameters.seed);
        this._is_initialized = false;
        return this;
    }

    /**
     * Set and get parameters
     * @param {string} [name = null] - Name of the parameter. If not given then returns all parameters as an Object.
     * @param {any} [value = null] - Value of the parameter to set. If <code>name</code> is set and <code>value</code> is not given, returns the value of the respective parameter.
     * @returns {DR|any|object}
     * On setting a parameter, this function returns the DR object.
     * If <code>name</code> is set and <code>value == null</code> then return actual parameter value.
     * If <code>name</code> is not given, then returns all parameters as an Object.
     *
     * @example
     * '''
     * const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
     * DR.parameter("d"); // returns 3,
     * DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
     * '''
     */
    parameter(name = null, value = null) {
        if (name === null) {
            return Object.assign({}, this._parameters);
        }
        if (!this._parameters.hasOwnProperty(name)) {
            throw new Error(`${name} is not a valid parameter!`);
        }
        if (value !== null) {
            this._parameters[name] = value;
            this._is_initialized = false;
            return this;
        } else {
            return this._parameters[name];
        }
    }

    para(name = null, value = null) {
        return this.parameter(name, value);
    }

    p(name = null, value = null) {
        return this.parameter(name, value);
    }

    /**
     * Computes the projection.
     * @returns {Matrix} the projection.
     */
    transform() {
        this.check_init();
        return this.projection;
    }

    /**
     * Computes the projection.
     * @yields {Matrix|number[][]} the intermediate steps of the projection.
     */
    *generator() {
        yield this.transform();
    }

    /**
     * If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.
     * @returns {DR}
     */
    check_init() {
        if (!this._is_initialized && typeof this.init === "function") {
            this.init();
            this._is_initialized = true;
        }
        return this;
    }

    /**
     * @returns {Matrix|number[][]} the projection in the type of input <code>X</code>.
     */
    get projection() {
        if (this.hasOwnProperty("Y")) {
            this.check_init();
            return this._type === "matrix" ? this.Y : this.Y.to2dArray;
        } else {
            throw new Error("The dataset is not transformed yet!");
        }
    }

    /**
     * Computes the projection.
     * @param  {...unknown} args - Arguments the transform method of the respective DR method takes.
     * @returns {Promise<Matrix|number[][]>} the dimensionality reduced dataset.
     */
    async transform_async(...args) {
        return this.transform(...args);
    }

    /**
     * Computes the projection.
     * @static
     * @param  {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Matrix|number[][]} the dimensionality reduced dataset.
     */
    static transform(...args) {
        let dr = new this(...args);
        return dr.transform();
    }

    /**
     * Computes the projection.
     * @static
     * @param  {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Promise} a promise yielding the dimensionality reduced dataset.
     */
    static async transform_async(...args) {
        return this.transform(...args);
    }

    /**
     * Computes the projection.
     * @static
     * @param  {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Generator} a generator yielding the intermediate steps of the dimensionality reduction method.
     */
    static *generator(...args) {
        const dr = new this(...args);
        const generator = dr.generator();
        for (const result of generator) {
            yield result;
        }
    }
}
