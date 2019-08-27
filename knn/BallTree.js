import { euclidean } from "../metrics/index";
import { Heap } from "../datastructure/index";
/**
 * @memberof module:knn
 */
export class BallTree {
    /**
     * Generates a BallTree with given {@link elements}.
     * @param {Array} elements - Elements which should be added to the BallTree
     * @param {function} metric metric to use: (a, b) => distance
     * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
     * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
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

    add(elements) {
        elements = elements.map((element, index) => {
            return {index: index, element: element}
        })
        this._root = this._construct(elements)
        return this;
    }

    _construct(elements) {
        if (elements.length === 1) {
            return new this._Leaf(elements);
        } else {
            let c = this._greatest_spread(elements);
            let sorted_elements = elements.sort((a, b) => a.element[c] - b.element[c])
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
                B = new this._Leaf(elements)
            }
            //if (this._root === null) this._root = B;
            return B;
        } /*else {
            return null;
        }*/
    }

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
        }, start)
        spread = spread.map(d => d[1] - d[0]);
        
        let c = 0;
        for (let i = 0; i < d; ++i) {
            c = spread[i] > spread[c] ? i : c;
        }
        return c
    }

    search(t, k=5) {
        return this._search(t, k, new Heap(null, d => this._metric(d.element, t), "max"), this._root);
    }

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
