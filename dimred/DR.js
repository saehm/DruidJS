import { euclidean } from "../metrics/index.js";
import { Matrix } from "../matrix/index.js";
import { Randomizer } from "../util/index.js";

/**
 * @class
 * @alias DR
 * @borrows DR#parameter as DR#para
 * @borrows DR#parameter as DR#p
 */
export class DR {
    //static parameter_list = [];
    get parameter_list() {
        return this._parameter_list;
    }

    set parameter_list(list) {
        this._parameter_list = list;
        return this;
    }
    /**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias DR
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.
     * @param {seed} [seed = 1212] - the seed value for the random number generator.
     * @returns {DR}
     */
    constructor(X, d = 2, metric = euclidean, seed = 1212) {
        if (Array.isArray(X)) {
            this._type = "array";
            this.X = Matrix.from(X);
        } else if (X instanceof Matrix) {
            this._type = "matrix";
            this.X = X;
        } else {
            throw new Error("no valid type for X");
        }
        [this._N, this._D] = this.X.shape;
        this._d = d;
        this._metric = metric;
        this._seed = seed;
        this._randomizer = new Randomizer(seed);
        this._is_initialized = false;
        return this;
    }

    /**
     * Set and get parameters
     * @param {String} name - name of the parameter.
     * @param {Number} [value = null] - value of the parameter to set, if <code>value == null</code> then return actual parameter value.
     * @memberof DR
     */
    parameter(name, value = null) {
        if (!this.parameter_list.includes(name)) {
            throw new Error(`${name} is not a valid parameter!`);
        }
        if (value) {
            this[`_${name}`] = value;
            this._is_initialized = false;
            return this;
        } else {
            return this[`_${name}`];
        }
    }

    para(name, value = null) {
        return this.parameter(name, value);
    }

    p(name, value = null) {
        return this.parameter(name, value);
    }

    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    transform() {
        this.check_init();
        return this.projection;
    }

    *generator() {
        return this.transform();
    }

    check_init() {
        if (!this._is_initialized && typeof this.init === "function") {
            this.init();
            this._is_initialized = true;
        }
    }

    /**
     * @returns {Matrix} Returns the projection.
     */
    get projection() {
        return this._type === "matrix" ? this.Y : this.Y.to2dArray;
    }

    async transform_async(...args) {
        const dr = new this(...args);
        return dr.transform();
    }

    static transform(...args) {
        let dr = new this(...args);
        return dr.transform();
    }

    static async transform_async(...args) {
        return this.transform(...args);
    }

    static *generator(...args) {
        const dr = new this(...args);
        const gen = dr.generator();
        for (const res of gen) {
            yield res;
        }
    }
}
