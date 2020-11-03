import { euclidean } from "../metrics/index";
import { Matrix } from "../matrix/index";
import { Randomizer } from "../util/randomizer";

/**
 * @class
 * @alias DR
 */
export class DR{
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
     * @param {seed} [seed=1987] - the seed value for the random number generator.
     * @returns {DR}
     */
    constructor(X, d=2, metric=euclidean, seed=1212) {
        if (Array.isArray(X)) {
            this._type = "array";
            this.X = Matrix.from(X);
        } else if (X instanceof Matrix) {
            this._type = "matrix";
            this.X = X;
        } else {
            throw "no valid type for X";
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
     * @param {Number} [value = null] - value of the parameter to set, if null then return actual parameter value.
     */
    parameter(name, value=null) {
        if (this.parameter_list.findIndex(parameter => parameter === name) === -1) {
            throw `${name} is not a valid parameter!`;
        } 
        if (value) {
            this[`_${name}`] = value;
            return this; 
        } else {
            return this[`_${name}`];
        }
    }

    /**
     * Alias for 'parameter'.
     * @param {String} name 
     * @param {Number} value 
     */
    para(name, value=null) {
        return this.parameter(name, value);
    }

    /**
     * Alias for 'parameter'.
     * @param {String} name 
     * @param {Number} value 
     */
    p(name, value=null) {
        return this.parameter(name, value);
    }

    /**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */
    transform() {
        this.check_init();
        return this.Y;
    }

    generator() {
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

    async transform_async() {
        return this.transform();
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