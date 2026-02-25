import { Matrix } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";

/** @import {InputType} from "../index.js" */

/**
 * @abstract
 * @template {InputType} T
 * @template {{ seed?: number }} Para
 *
 * Base class for all Dimensionality Reduction (DR) algorithms.
 *
 * Provides a common interface for parameters management, data initialization,
 * and transformation (both synchronous and asynchronous).
 *
 * @class
 */
export class DR {
    /** @type {number} */
    _D;
    /** @type {number} */
    _N;
    /** @type {Randomizer} */
    _randomizer;
    /** @type {boolean} */
    _is_initialized;

    /**
     * Takes the default parameters and seals them, remembers the type of input `X`, and initializes the random number
     * generator.
     *
     * @param {T} X - The high-dimensional data.
     * @param {Para} default_parameters - Object containing default parameterization of the DR method.
     * @param {Partial<Para>} parameters - Object containing parameterization of the DR method to override defaults.
     */
    constructor(X, default_parameters, parameters = {}) {
        /** @type {T} */
        this.__input = X;

        /** @type {Para} */
        this._parameters = /** @type {Para} */ Object.seal({
            ...default_parameters,
            ...parameters,
        });
        /** @type {"array" | "matrix" | "typed"} */
        this._type;
        /** @type {Matrix} */
        this.X;
        /** @type {Matrix} */
        this.Y;

        if (Array.isArray(X)) {
            if (X[0] instanceof Float64Array) {
                this._type = "typed";
            } else {
                this._type = "array";
            }
            this.X = Matrix.from(X);
        } else if (X instanceof Matrix) {
            this._type = "matrix";
            this.X = X;
        } else {
            throw new Error("No valid type for X!");
        }
        const [N, D] = this.X.shape;
        this._N = N;
        this._D = D;
        this._randomizer = new Randomizer(this._parameters.seed);
        this._is_initialized = false;
    }

    /**
     * Get all Parameters.
     * @overload
     * @returns {Para}
     */
    /**
     * Get value of given parameter.
     * @template {keyof Para} K
     * @overload
     * @param {K} name - Name of the parameter.
     * @returns {Para[K]}
     */
    /**
     * Set value of given parameter.
     * @template {keyof Para} K
     * @overload
     * @param {K} name - Name of the parameter.
     * @param {Para[K]} value - Value of the parameter to set.
     * @returns {this}
     */
    /**
     * @param {keyof Para} [name] - Name of the parameter. If null, returns all parameters as an Object.
     * @param {Para[keyof Para]} [value] - Value of the parameter to set. If name is set and value is not given, returns the
     *   current value.
     * @returns {Para | Para[keyof Para] | this} On setting a parameter, returns the DR object. If name is set and value is not
     *   given, returns the parameter value. If name is null, returns all parameters. On setting a parameter, this
     *   function returns the DR object. If `name` is set and `value == null` then return actual parameter value. If
     *   `name` is not given, then returns all parameters as an Object.
     * @example
     * ```js
     * const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for `d = 3`.
     * DR.parameter("d"); // returns 3
     * DR.parameter("d", 2); // sets parameter `d` to 2 and returns `DR`.
     * ```
     *
     */
    parameter(name, value) {
        if (name === undefined && value === undefined) {
            return Object.assign({}, this._parameters);
        }
        if (name && !Object.hasOwn(this._parameters, name)) {
            throw new Error(`${String(name)} is not a valid parameter!`);
        }
        if (name && value !== undefined) {
            this._parameters[name] = value;
            this._is_initialized = false;
            return this;
        } else if (name) {
            return this._parameters[name];
        }
        throw new Error("Should not happen!");
    }

    /**
     * Computes the projection.
     *
     * @abstract
     * @param {...unknown} args
     * @returns {T} The projection.
     */
    transform(...args) {
        args;
        this.check_init();
        return this.projection;
    }

    /**
     * Computes the projection.
     *
     * @template {InputType} T
     * @template {{ seed?: number }} Para
     * @param {T} X
     * @param {Para} parameters
     * @param {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {T} The dimensionality reduced dataset.
     */
    static transform(X, parameters, ...args) {
        args;
        const dr = new DR(X, parameters, parameters);
        return /** @type {T} */ (dr.transform());
    }

    /**
     * Computes the projection.
     *
     * @abstract
     * @param {...unknown} args
     * @returns {Generator<T, T, void>} The intermediate steps of the projection.
     */
    *generator(...args) {
        const R = this.transform(...args);
        yield R;
        return R;
    }

    /**
     * Computes the projection.
     *
     * @template {{ seed?: number }} Para
     * @param {InputType} X
     * @param {Para} parameters
     * @param {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Generator<InputType, InputType, void>} A generator yielding the intermediate steps of the dimensionality
     *   reduction method.
     */
    static *generator(X, parameters, ...args) {
        const dr = new DR(X, parameters, parameters);
        const generator = dr.generator(...args);
        let result;
        do {
            result = generator.next();
            yield result.value;
        } while (!result.done);

        return result.value;
    }

    /**
     * @abstract
     * @param {...unknown} args
     */
    init(...args) {
        args;
    }

    /**
     * If the respective DR method has an `init` function, call it before `transform`.
     *
     * @returns {DR<T, Para>}
     */
    check_init() {
        if (!this._is_initialized && typeof this.init === "function") {
            this.init();
            this._is_initialized = true;
        }
        return this;
    }

    /** @returns {T} The projection in the type of input `X`. */
    get projection() {
        if (Object.hasOwn(this, "Y")) {
            this.check_init();
            //return this._type === "matrix" ? this.Y : this.Y.to2dArray();
            if (this._type === "matrix") {
                return /** @type {T} */ (/** @type {any} */ (this.Y));
            } else if (this._type === "typed") {
                return /** @type {T} */ (/** @type {any} */ (this.Y.to2dArray()));
            } else {
                return /** @type {T} */ (/** @type {any} */ (this.Y.asArray()));
            }
        } else {
            throw new Error("The dataset is not transformed yet!");
        }
    }

    /**
     * Computes the projection.
     *
     * @param {...unknown} args - Arguments the transform method of the respective DR method takes.
     * @returns {Promise<T>} The dimensionality reduced dataset.
     */
    async transform_async(...args) {
        return this.transform(...args);
    }

    /**
     * Computes the projection.
     *
     * @template {{ seed?: number }} Para
     * @param {InputType} X
     * @param {Para} parameters
     * @param {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Promise<X>} A promise yielding the dimensionality reduced dataset.
     */
    static async transform_async(X, parameters, ...args) {
        return DR.transform(X, parameters, ...args);
    }
}
