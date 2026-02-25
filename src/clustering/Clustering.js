import { Matrix } from "../matrix/index.js";

/** @import {InputType} from "../index.js" */

/**
 * Base class for all clustering algorithms.
 * @template Para
 */
export class Clustering {
    /** @type {InputType} */
    _points;
    /** @type {Para} */
    _parameters;
    /** @type {Matrix} */
    _matrix;
    /** @type {number} */
    _N;
    /** @type {number} */
    _D;

    /**
     * Compute the respective Clustering with given parameters
     * @param {InputType} points
     * @param {Para} parameters
     */
    constructor(points, parameters) {
        this._points = points;
        this._parameters = parameters;

        this._matrix = points instanceof Matrix ? points : Matrix.from(points);
        const [N, D] = this._matrix.shape;
        this._N = N;
        this._D = D;
    }

    /**
     * @abstract
     * @param {...unknown} args
     * @returns {number[][]} An array with the indices of the clusters.
     */
    get_clusters(...args) {
        args;
        throw new Error("The function get_clusters must be implemented!");
    }

    /**
     * @abstract
     * @param {...unknown} args
     * @returns {number[]} An array with the clusters id's for each point.
     */
    get_cluster_list(...args) {
        args;
        throw new Error("The function get_cluster_list must be implemented!");
    }
}
