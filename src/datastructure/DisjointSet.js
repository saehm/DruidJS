/**
 * @template T
 * @typedef {Object} DisjointSetPayload
 * @property {T} parent
 * @property {Set<T>} children
 * @property {number} size
 */

/**
 * @template T
 * @class
 * @category Data Structures
 * @see {@link https://en.wikipedia.org/wiki/Disjoint-set_data_structure}
 */
export class DisjointSet {
    /**
     * @param {T[]?} elements
     */
    constructor(elements = null) {
        /**
         * @private
         * @type {Map<T, DisjointSetPayload<T>>}
         */
        this._list = new Map();
        if (elements) {
            for (const e of elements) {
                this.make_set(e);
            }
        }
    }

    /**
     * @private
     * @param {T} x
     * @returns {DisjointSet<T>}
     */
    make_set(x) {
        const list = this._list;
        if (!list.has(x)) {
            list.set(x, { parent: x, children: new Set([x]), size: 1 });
        }
        return this;
    }

    /**
     * @param {T} x
     * @returns
     */
    find(x) {
        const list = this._list;
        const disjoint_set = list.get(x);
        if (disjoint_set) {
            if (disjoint_set.parent !== x) {
                disjoint_set.children.add(x);
                const new_parent = this.find(disjoint_set.parent);
                if (!new_parent) throw new Error("should not happen!");
                disjoint_set.parent = new_parent;
                return disjoint_set.parent;
            } else {
                return x;
            }
        } else {
            return null;
        }
    }

    /**
     * @param {T} x
     * @param {T} y
     * @returns
     */
    union(x, y) {
        let node_x = this.find(x);
        let node_y = this.find(y);

        if (!node_x || !node_y) throw new Error("x or y not found!");

        let disjoint_set_x = this._list.get(node_x);
        let disjoint_set_y = this._list.get(node_y);

        if (!disjoint_set_x || !disjoint_set_y) throw new Error("should not happen!");

        if (node_x === node_y) return this;
        if (disjoint_set_x.size < disjoint_set_y.size) {
            [node_x, node_y] = [node_y, node_x];
            [disjoint_set_x, disjoint_set_y] = [disjoint_set_y, disjoint_set_x];
        }

        disjoint_set_y.parent = node_x;
        // keep track of children
        disjoint_set_y.children.forEach(disjoint_set_x.children.add, disjoint_set_x.children);
        disjoint_set_x.size += disjoint_set_y.size;

        return this;
    }

    /** @param {T} x */
    get_children(x) {
        const node = this._list.get(x);
        if (node) {
            return node.children;
        } else {
            return null;
        }
    }
}
