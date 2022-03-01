/**
 * @class
 * @alias KNN
 */
export class KNN {
    /**
     * Generates a KNN list with given {@link elements}.
     * @constructor
     * @memberof module:knn
     * @alias KNN
     * @param {Array=} elements - Elements which should be added to the KNN list
     * @param {Function|"precomputed"} [metric = euclidean] metric is either precomputed or a function to use: (a, b) => distance
     * @returns {KNN}
     */
    constructor(elements?: any[] | undefined, metric?: Function | "precomputed");
    _metric: Function | "precomputed";
    _elements: Matrix;
    _D: Matrix;
    KNN: Heap[];
    /**
     *
     * @param {Array|Number} t - query element or index.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @returns {Heap} - Heap consists of the {@link k} nearest neighbors.
     */
    search(t: any[] | number, k?: number): Heap;
}
import { Matrix } from "../matrix/index.js";
import { Heap } from "../datastructure/index.js";
//# sourceMappingURL=KNN.d.ts.map