import { euclidean } from '../metrics/index.js';
import { Heap } from '../datastructure/index.js';
import { Randomizer } from '../util/index.js';

/**
 * @class
 * @alias HNSW
 */
export class HNSW {
    /**
     * Hierarchical navigable small world graph. Efficient and robust approximate nearest neighbor search.
     * @constructor
     * @memberof module:knn
     * @alias HNSW
     * @param {Function} [metric = euclidean] - metric to use: (a, b) => distance.
     * @param {Boolean} [heuristic = true] - use heuristics or naive selection.
     * @param {Number} [m = 5] - max number of connections.
     * @param {Number} [ef = 200] - size of candidate list.
     * @param {Number} [m0 = 2 * m] - max number of connections for ground layer.
     * @param {Number} [mL = 1 / Math.log2(m)] - normalization factor for level generation.
     * @param {Number} [seed = 1987] - seed for random number generator.
     * @see {@link https://arxiv.org/abs/1603.09320}
     * @see {@link https://arxiv.org/pdf/1904.02077}
     */
    constructor(metric = euclidean, heuristic = true, m = 5, ef = 200, m0 = null, mL = null, seed = 1987) {
        this._metric = metric;
        this._select = heuristic ? this._select_heuristic : this._select_simple;
        this._m = m;
        this._ef = ef;
        this._m0 = m0 || 2 * m;
        this._graph = new Map();
        this._ep = null;
        this._L = null;
        this._mL = mL || 1 / Math.log2(m);
        this._randomizer = new Randomizer(seed);
    }

    addOne(element) {
        this.add([element])
    }

    /**
     * 
     * @param {Array<*>} elements - new elements.
     * @returns {HNSW}
     */
    add(elements) {
        const m = this._m;
        const ef = this._ef;
        const m0 = this._m0;
        const mL = this._mL;
        const randomizer = this._randomizer;
        let graph = this._graph;
        for (const element of elements) {
            let ep = this._ep ? this._ep.slice() : null;
            let W = [];
            const L = this._L;
            const rand = Math.min(randomizer.random + 1e-8, 1);
            let l = Math.floor(-Math.log(rand * mL))
            const min_L_l = Math.min(L, l);
            if (L) {
                for (let l_c = graph.size - 1; l_c > min_L_l; --l_c) {
                    ep = this._search_layer(element, ep, 1, l_c);
                }
                for (let l_c = min_L_l; l_c >= 0; --l_c) {
                    const layer_c = graph.get(l_c);//[l_c];
                    layer_c.points.push(element)
                    W = this._search_layer(element, ep, ef, l_c);
                    const neighbors = l_c > 3 ? this._select(element, W, m, l_c) : this._select_simple(element, W, m);
                    for (const p of neighbors) {
                        if (p !== element) {
                            const edges_p = layer_c.edges.get(p);
                            if (!edges_p) {
                                layer_c.edges.set(p, [element]) ;
                            }else {
                                edges_p.push(element);
                            }
                            const edges_e = layer_c.edges.get(element);
                            if (!edges_e) {
                                layer_c.edges.set(element, [p]) ;
                            } else {
                                edges_e.push(p);
                            }
                        }
                    }
                    const max = (l_c === 0 ? m0 : m);
                    for (const e of neighbors) {
                        const e_conn = layer_c.edges.get(e);
                        if (e_conn.length > max) {
                            const neighborhood = this._select(e, e_conn, max, l_c);
                            layer_c.edges.delete(e);
                            layer_c.edges.set(e, neighborhood);
                        }
                    }
                    ep = W;
                }
            }
            let N = graph.size;
            if (N < l || l > L) {
                for (let i = N; i <= l; ++i) {
                    graph.set(i, {
                        "l_c": i, 
                        "points": [element], 
                        "edges": new Map(),
                    });
                }
                this._ep = [element];
                this._L = l;
            }
        }
        return this;
    }

    /**
     * @private
     * @param {*} q - base element.
     * @param {Array} candidates - candidate elements.
     * @param {Number} M - number of neighbors to return.
     * @param {Number} l_c - layer number.
     * @param {Boolean} [extend_candidates = true] - flag indicating wheter or not to extend candidate list.
     * @param {Boolean} [keep_pruned_connections = true] - flag indicating wheter or not to add discarded elements.
     * @returns M elements selected by the heuristic.
     */
    _select_heuristic(q, candidates, M, l_c, extend_candidates = true, keep_pruned_connections = true) {
        if (l_c > this._graph.size - 1) return candidates
        const metric = this._metric;
        const randomizer = this._randomizer;
        const layer = this._graph.get(l_c);
        let R = [];
        let W_set = new Set(candidates);
        if (extend_candidates) {
            for (const c of candidates) {
                const edges = layer.edges.get(c);
                if (!edges) break;
                for (const c_adj of edges) {
                    W_set.add(c_adj)
                }
            }
        }
        let W = new Heap(W_set, d => metric(d, q), "min")
        let W_d = new Heap(null, d => metric(d, q), "min");
        while (!W.empty && R.length < M) {
            let e = W.pop()
            let random_r = randomizer.random_int % R.length;
            if (R.length === 0 || e.value < metric(R[random_r], q)) {
                R.push(e.element);
            } else {
                W_d.push(e.element)
            }
        }
        if (keep_pruned_connections) {
            while (!W_d.empty && R.length < M) {
                R.push(W_d.pop().element)
            }
        }
        return R
    }

    /**
     * @private
     * @param {*} q - base element.
     * @param {Array} C - candidate elements.
     * @param {Number} M - number of neighbors to return.
     * @returns {Array} M nearest elements from C to q.
     */
    _select_simple(q, C, M) {
        const metric = this._metric;
        let res = C.sort((a,b) => metric(a, q) - metric(b, q)).slice(0,M);
        return res
    }

    /**
     * @private
     * @param {*} q - query element.
     * @param {Array} ep - enter points.
     * @param {Number} ef - number of nearest to {@link q} elements to return.
     * @param {Number} l_c - layer number.
     * @returns {Array} ef closest neighbors to q.
     */
    _search_layer(q, ep, ef, l_c) {
        const metric = this._metric;
        const layer = this._graph.get(l_c)//.find(l => l.l_c === l_c);//[l_c];
        if (layer.edges.size === 0) return ep;
        let v = new Set(ep);
        let C = new Heap(v, d => metric(d, q), "min");
        let W = new Heap(v, d => metric(d, q), "max");
        while (!C.empty) {
            const c = C.pop();
            let f = W.first;
            if (c.value > f.value) {
                break;
            }
            const edges = layer.edges.get(c.element);
            if (!edges) break;
            for (const e of edges) {
                if (!v.has(e)) {
                    v.add(e);
                    f = W.first;
                    if (metric(e, q) < metric(f.element, q) || W.length < ef) {
                        C.push(e);
                        W.push(e);
                        if (W.length > ef) {
                            W.pop();
                        }
                    }
                }
            }
        }
        return W.data();
    }

    /**
     * 
     * @param {*} q - query element.
     * @param {*} K - number of nearest neighbors to return.
     * @param {*} ef - size of the dynamic cnadidate list.
     * @returns {Array} K nearest elements to q.
     */
    search(q, K, ef = 1) {
        let ep = this._ep.slice();
        let L = this._L;
        for (let l_c = L; l_c > 0; --l_c) {
            ep = this._search_layer(q, ep, ef, l_c);
        }
        ep = this._search_layer(q, ep, K, 0);
        return ep;
    }

    /**
     * Iterator for searching the HNSW graphs
     * @param {*} q - query element.
     * @param {*} K - number of nearest neighbors to return.
     * @param {*} ef - size of the dynamic cnadidate list.
     * @yields {Array} K nearest elements to q.
     */
    * search_iter(q, K, ef = 1) {
        let ep = this._ep.slice();
        let L = this._L;
        yield{"l_c": L, "ep": [q]}
        for (let l_c = L; l_c > 0; --l_c) {
            yield {"l_c": l_c, "ep": ep}
            ep = this._search_layer(q, ep, ef, l_c);
            yield {"l_c": l_c, "ep": ep}
        }
        yield {"l_c": 0, "ep": ep}
        ep = this._search_layer(q, ep, K, 0);
        yield {"l_c": 0, "ep": ep}
    }
}