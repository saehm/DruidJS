/**
 * @class
 * @alias BallTree
 */
export class BallTree {
    /**
     * Generates a BallTree with given {@link elements}.
     * @constructor
     * @memberof module:knn
     * @alias BallTree
     * @param {Array=} elements - Elements which should be added to the BallTree
     * @param {Function} [metric = euclidean] metric to use: (a, b) => distance
     * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
     * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
     * @returns {BallTree}
     */
    constructor(elements?: any[] | undefined, metric?: Function);
    _Node: {
        new (pivot: any, child1?: any, child2?: any, radius?: any): {
            pivot: any;
            child1: any;
            child2: any;
            radius: any;
        };
    };
    _Leaf: {
        new (points: any): {
            points: any;
        };
    };
    _metric: Function;
    /**
     *
     * @param {Array<*>} elements - new elements.
     * @returns {BallTree}
     */
    add(elements: Array<any>): BallTree;
    _root: Node;
    /**
     * @private
     * @param {Array<*>} elements
     * @returns {Node} root of balltree.
     */
    private _construct;
    /**
     * @private
     * @param {Node} B
     * @returns {Number}
     */
    private _greatest_spread;
    /**
     *
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @returns {Heap} - Heap consists of the {@link k} nearest neighbors.
     */
    search(t: any, k?: number): Heap;
    /**
     * @private
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @param {Heap} Q - Heap consists of the currently found {@link k} nearest neighbors.
     * @param {Node|Leaf} B
     */
    private _search;
}
import { Heap } from "../datastructure/index.js";
//# sourceMappingURL=BallTree.d.ts.map