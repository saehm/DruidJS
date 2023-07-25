import { euclidean } from "../metrics/index.js";
import { Heap } from "../datastructure/index.js";
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
    constructor(elements = null, metric = euclidean) {
        this._Node = class {
            constructor(pivot, child1=null, child2=null, radius=null) {
                this.pivot = pivot;
                this.child1 = child1;
                this.child2 = child2;
                this.radius = radius;
            }
        }
        this._Leaf = class {
            constructor(points) {
                this.points = points;
            }
        }
        this._metric = metric;
        if (elements) {
            this.add(elements);
        }
        return this;
    }

    /**
     * 
     * @param {Array<*>} elements - new elements.
     * @returns {BallTree}
     */
    add(elements) {
        elements = elements.map((element, index) => {
            return {index: index, element: element}
        })
        this._root = this._construct(elements);
        return this;
    }

    /**
     * @private
     * @param {Array<*>} elements 
     * @returns {Node} root of balltree.
     */
    _construct(elements) {
        if (elements.length === 1) {
            return new this._Leaf(elements);
        } else {
            let c = this._greatest_spread(elements);
            let sorted_elements = elements.sort((a, b) => a.element[c] - b.element[c]);
            let n = sorted_elements.length;
            let p_index = Math.floor(n / 2);
            let p = elements[p_index];
            let L = sorted_elements.slice(0, p_index);
            let R = sorted_elements.slice(p_index, n);
            let radius = Math.max(...elements.map(d => this._metric(p.element, d.element)));
            let B
            if (L.length > 0 && R.length > 0) {         
                B = new this._Node(p, this._construct(L), this._construct(R), radius);
            } else {
                B = new this._Leaf(elements);
            }
            return B;
        }
    }

    /**
     * @private
     * @param {Node} B 
     * @returns {Number}
     */
    _greatest_spread(B) {
        let d = B[0].element.length;
        let start = new Array(d);

        for (let i = 0; i < d; ++i) {
            start[i] = [Infinity, -Infinity];
        }

        let spread = B.reduce((acc, current) => {
            for (let i = 0; i < d; ++i) {
                acc[i][0] = Math.min(acc[i][0], current.element[i]);
                acc[i][1] = Math.max(acc[i][1], current.element[i]);
            }
            return acc;
        }, start);
        spread = spread.map(d => d[1] - d[0]);
        
        let c = 0;
        for (let i = 0; i < d; ++i) {
            c = spread[i] > spread[c] ? i : c;
        }
        return c;
    }

    /**
     * 
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @returns {Heap} - Heap consists of the {@link k} nearest neighbors.
     */
    search(t, k = 5) {
        return this._search(t, k, new Heap(null, d => this._metric(d.element, t), "max"), this._root);
    }

    /**
     * @private
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @param {Heap} Q - Heap consists of the currently found {@link k} nearest neighbors.
     * @param {Node|Leaf} B 
     */
    _search(t, k, Q, B) {
        // B is Node
        if (Q.length >= k && B.pivot && B.radius && this._metric(t, B.pivot.element) - B.radius >= Q.first.value) {
            return Q;
        } 
        if (B.child1) this._search(t, k, Q, B.child1);
        if (B.child2) this._search(t, k, Q, B.child2);
        
        // B is leaf
        if (B.points) {
            for (let i = 0, n = B.points.length; i < n; ++i) {
                let p = B.points[i];
                if (k > Q.length) {
                    Q.push(p);
                } else {
                    Q.push(p);
                    Q.pop();
                }
            }
        }
        return Q;
    }
}