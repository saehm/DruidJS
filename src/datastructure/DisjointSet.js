/**
 * @class
 * @alias DisjointSet
 * @see {@link https://en.wikipedia.org/wiki/Disjoint-set_data_structure}
 */
export class DisjointSet {
    /**
     * @constructor
     * @alias DisjointSet
     * @memberof module:datastructure
     * @param {Array=} elements 
     * @returns {DisjointSet}
     */
    constructor(elements = null) {
        this._list = new Set();
        if (elements) {
            for (const e of elements) {
                this.make_set(e);
            }
        }
        return this;
    }

    make_set(x) {
        const list = this._list;
        if (!list.has(x)) {
            list.add(x);
            x.__disjoint_set = {};
            x.__disjoint_set.parent = x;
            x.__disjoint_set.children = new Set([x]);
            x.__disjoint_set.size = 1;
        }
        return this;
    }

    find(x) {
        const list = this._list;
        if (list.has(x)) {
            if (x.__disjoint_set.parent !== x) {
                x.__disjoint_set.children.add(...x);
                x.__disjoint_set.parent = this.find(x.__disjoint_set.parent);
                return x.__disjoint_set.parent;
            } else {
                return x;
            }
        } else {
            return null;
        }
    }

    union(x, y) {
        let node_x = this.find(x);
        let node_y = this.find(y);

        if (node_x === node_y) return this;
        if (node_x.__disjoint_set.size < node_y.__disjoint_set.size) [node_x, node_y] = [node_y, node_x];

        node_y.__disjoint_set.parent = node_x;
        // keep track of children?
        node_y.__disjoint_set.children.forEach(node_x.__disjoint_set.children.add, node_x.__disjoint_set.children);
        node_x.__disjoint_set.size += node_y.__disjoint_set.size;

        return this;
    }
}