// https://null.org v0.0.1 Copyright 2019 abc
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(global = global || self, factory(global.visci = global.visci || {}));
}(this, function (exports) { 'use strict';

function zeros(n=0, m=0) {
    if (n < 0 || m < 0) {
        return undefined
    } else if (n === 0 && m === 0) {
        return 0
    } else if (n === 0 || m === 0) {
        let A = new Array(Math.max(n,m));
        for (let i = 0; i < m; ++i) {
            A[i] = 0;
        }
        return A
    } else {    
        let A = new Array(n);
        for (let i = 0; i < n; ++i) {
            A[i] = new Array(m);
            for (let j = 0; j < m; ++j) {
                A[i][j] = 0;
            }
        }
        return A
    }
}

function euclidean(a, b) {
    if (a.length != b.length) return undefined
    let n = a.length;
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += ((a[i] - b[i]) ** 2);
    }
    return Math.sqrt(sum)
}

function cosine(a, b) {
    if (a.length != b.length) return undefined;
    let n = a.length;
    let sum = 0;
    let sum_a = 0;
    let sum_b = 0;
    for (let i = 0; i < n; ++i) {
        sum += (a[i] * b[i]);
        sum_a += (a[i] * a[i]);
        sum_b += (b[i] * b[i]);
    }
    return sum / (Math.sqrt(sum_a) * Math.sqrt(sum_b));
}

const euclidean$1 = euclidean;

function dmatrix(A, metric = euclidean$1) {
    let distance = metric;
    if (distance === undefined) return undefined
    let n = A.length;
    let D = zeros(n,n);
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            D[i][j] = D[j][i] = distance(A[i], A[j]);
        }
    }
    return D
}

const euclidean$2 = euclidean;

function k_nearest_neighbors(A, k, distance_matrix = null, metric = euclidean$2) {
    let n = A.length;
    let D = distance_matrix || dmatrix(A, metric);
    for (let i = 0; i < n; ++i) {
        D[i] = D[i].map((d,j) => {
            return {
                i: i, j: j, distance: D[i][j]
            }
        }).sort((a, b) => a.distance - b.distance)
        .slice(1, k + 1);
    }
    return D
}

function linspace(start, end, number) {
    if (number === undefined) {
        number = Math.max(Math.round(end - start) + 1, 1);
    }
    if (number < 2) {
        return n === 1 ? [a] : [];
    }
    let result = new Array(number);
    number -= 1;
    for (let i = number; i >= 0; --i) {
        result[i] = (i * end + (number - i) * start) / number;
    }
    return result
}

class Heap {
    constructor(arr = null, accessor = (d) => d, comparator = "min") {
        this.root = null;
        this.accessor = accessor;

        if (comparator == "min") {
            this._comparator = (a, b) => a <= b;
        } else if (comparator == "max") {
            this._comparator = (a, b) => a >= b;
        } else {
            this._comparator = comparator;
        }

        if (arr && arr.length > 0) {
            let self = this;
            arr.forEach(d => self.push(d));
        }
    }

    push(element) {
        const value = this.accessor(element);
        const newNode = new Node(element, value);
        if (!this.root || this._comparator(value, this.root.value)) {
            newNode.next = this.root;
            this.root = newNode;
        } else {
            let pointer = this.root;
            while (pointer.next && !this._comparator(value, pointer.next.value)) {
                pointer = pointer.next;
            }
            newNode.next = pointer.next;
            pointer.next = newNode;
        }
        return this;
    }

    pop() {
        if (!this.root) {
            return null;
        }
        const root = this.root;
        this.root = this.root.next;
        return root;
    }

    get first() {
        return this.root;
    }

    *iterate() {
        let pointer = this.root;
        /*do {
            yield pointer.element;
        } while (pointer = pointer.next);*/
        while (pointer) {
            yield pointer.element;
            pointer = pointer.next;
        }
    }

    toArray() {
        let res = [];
        let pointer = this.root;
        while (pointer) {
            res.push(pointer.element);
            pointer = pointer.next;
        }
        return res;
    }

    get length() {
        let len = 0;
        let pointer = this.root;
        while (pointer) {
            len += 1;
            pointer = pointer.next;
        }
        return len;
    }

    get empty() {
        return this.root === null;
    }
}

class Node {
    constructor(element, value) {
        this.element = element;
        this.value = value;
        this.next = null;
    }
}

const euclidean$3 = euclidean;

class HNSW {
    /**
     * 
     * @param {*} metric metric to use: (a, b) => distance
     * @param {*} heuristic use heuristics or naive selection
     * @param {*} m max number of connections
     * @param {*} ef size of candidate list
     * @param {*} m0 max number of connections for ground layer 
     */
    constructor(metric = euclidean$3, heuristic = true, m = 5, ef = 200, m0 = null) {
        this._metric = metric;
        this._select = heuristic ? this._select_heuristic : this._select_simple;
        this._m = m;
        this._ef = ef;
        this._m0 = m0 || 2 * m;
        this._graph = [];
        this._ep = null;
        this._L = null;
        this._mL = 1 / Math.log2(m);
        this.search = this.search;
    }

    addOne(element) {
        this.add([element]);
    }

    add(...elements) {
        const m = this._m;
        const ef = this._ef;
        const m0 = this._m0;
        //const metric = this._metric;
        const mL = this._mL;
        let graph = this._graph;
        for (const element of elements) {
            let ep = this._ep;
            let W = [];
            let L = this._L;
            let l = Math.floor(-Math.log(Math.random() * mL)) + 1;
            let min_L_l = Math.min(L, l);
            if (L) {
                for (let l_c = graph.length - 1; l_c > min_L_l; --l_c) {
                    W = this._search_layer(element, ep, 1, l_c);
                    ep = W;
                }
                for (let l_c = min_L_l; l_c >= 0; --l_c) {
                    let layer_c = graph[l_c];
                    layer_c.points.push(element);
                    W = this._search_layer(element, ep, ef, l_c);
                    let neighbors = this._select(element, W, m, l_c);
                    neighbors.forEach(p => {
                        if (p !== element) {
                            //let distance = metric(p, element);
                            layer_c.edges.push({
                                idx1: p, 
                                idx2: element, 
                                ///distance: distance
                            });
                            layer_c.edges.push({
                                idx1: element, 
                                idx2: p, 
                                //distance: distance
                            });
                        }
                    });
                    let max = (l_c === 0 ? m0 : m);
                    for (let e of neighbors) {
                        let e_conn = layer_c.edges
                            .filter(edge => edge.idx1 === e)
                            .map(edge => edge.idx2);
                        if (e_conn.length > max) {
                            let neighborhood = this._select(e, e_conn, max, l_c);
                            layer_c.edges = layer_c.edges
                                .filter(edge => edge.idx1 !== e);
                            neighborhood.forEach(neighbor => {
                                if (e !== neighbor) {
                                    //let distance = metric(e, neighbor);
                                    layer_c.edges.push({
                                        idx1: e, 
                                        idx2: neighbor, 
                                        //distance: distance
                                    });
                                }
                            });
                        }
                    }
                    ep = W;
                }
            }
            if (graph.length < l || l > L) {
                for (let i = l, n = graph.length; i >= n; --i) {
                    let new_layer = {
                        l_c: i, 
                        points: [element], 
                        edges: []
                    };
                    graph.push(new_layer);
                    if (i === l) {
                        this._ep = new_layer.points;
                        this._L = l;
                    }
                }
                graph = graph.sort((a, b) => a.l_c - b.l_c);
            }
        }
        return this;
    }

    _select_heuristic(q, candidates, M, l_c, extend_candidates = true, keep_pruned_connections = true) {
        if (l_c > this._graph.length - 1) return candidates
        const metric = this._metric;
        const layer = this._graph[l_c];
        let R = [];
        let W_set = new Set();
        candidates.forEach(c => W_set.add(c));
        if (extend_candidates) {
            for (let c of candidates) {
                for (let {idx2: c_adj} of layer.edges.filter(edge => edge.idx1 === c)) {
                    W_set.add(c_adj);
                }
            }
        }
        let W = new Heap(Array.from(W_set), d => metric(d, q), "min");
        let W_d = new Heap(null, d => metric(d, q), "min");
        while (W.first && R.length < M) {
            let e = W.pop();
            let random_r = Math.floor(Math.random() * R.length);
            if (R.length === 0 || e.value < metric(R[random_r], q)) {
                R.push(e.element);
            } else {
                W_d.push(e.element);
            }
        }
        if (keep_pruned_connections) {
            while (W_d.first && R.length < M) {
                R.push(W_d.pop().element);
            }
        }
        return R
    }

    _select_simple(q, C, M) {
        const metric = this._metric;
        let res = C.sort((a,b) => metric(a, q) - metric(b, q)).slice(0,M);
        return res
    }

    _search_layer(q, ep, ef, l_c) {
        const metric = this._metric;
        const layer = this._graph[l_c];
        let v = new Set();
        ep.forEach(p => v.add(p));
        let C = new Heap(ep, d => metric(d, q), "min");
        let W = new Heap(ep, d => metric(d, q), "max");
        while (C.length > 0) {
            let c = C.pop().element;
            let f = W.first;
            if (c.value > f.value) {
                break;
            }
            for (let {idx2: e} of layer.edges.filter(e => e.idx1 === c)) {
                if (!v.has(e)) {
                    v.add(e);
                    f = W.first.element;
                    if (metric(e, q) < metric(f, q) || W.length < ef) {
                        C.push(e);
                        W.push(e);
                        if (W.length > ef) {
                            W.pop();
                        }
                    }
                }
            }
        }
        let res =  W.toArray().reverse();
        return res;
    }

    search(q, K, ef = null) {
        ef = ef || this._ef;
        let ep = this._ep;
        let L = this._L;
        for (let l_c = L; l_c > 0; --l_c) {
            ep = this._search_layer(q, ep, 1, l_c);
        }
        ep = this._search_layer(q, ep, ef, 0);
        return ep.slice(0,K);
    }
}

exports.k_nearest_neighbors = k_nearest_neighbors;
exports.distance_matrix = dmatrix;
exports.zeros = zeros;
exports.linspace = linspace;
exports.HNSW = HNSW;
exports.Heap = Heap;
exports.euclidean = euclidean;
exports.cosine = cosine;

Object.defineProperty(exports, '__esModule', { value: true });

}));
