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
    constructor(elements?: any[] | undefined);
    _list: Set<any>;
    make_set(x: any): DisjointSet;
    find(x: any): any;
    union(x: any, y: any): DisjointSet;
}
//# sourceMappingURL=DisjointSet.d.ts.map