import { Heap } from "../datastructure/index.js";
import { euclidean } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { KNN } from "./KNN.js";

/** @import {ParametersNNDescent} from "./index.js" */
/**
 *
 * @template {number[] | Float64Array} T
 * @typedef {Object} NNDescentElement
 * @property {T} value
 * @property {number} index
 * @property {boolean} flag
 */

/**
 * @template {number[] | Float64Array} T
 * @typedef {Object} NNDescentNeighbor
 * @property {T} value
 * @property {number} index
 * @property {number} distance
 * @property {boolean} [flag]
 */

/**
 * NN-Descent
 *
 * An efficient graph-based approximate nearest neighbor search algorithm.
 * It works by iteratively improving a neighbor graph using the fact that
 * "neighbors of neighbors are likely to be neighbors".
 *
 * @class
 * @category KNN
 * @template {number[] | Float64Array} T
 * @extends KNN<T, ParametersNNDescent>
 * @see {@link http://www.cs.princeton.edu/cass/papers/www11.pdf|NN-Descent Paper}
 */
export class NNDescent extends KNN {
    /**
     * @private
     * @type {KNNHeap<T>[]}
     */
    _B = [];
    /**
     * @private
     * @type {NNDescentNeighbor<T>[][]}
     */
    nn = [];

    /**
     * @param {T[]} elements - Called V in paper.
     * @param {Partial<ParametersNNDescent>} [parameters={}] - Anything left out falls back to the
     *   documented default.
     * @see {@link http://www.cs.princeton.edu/cass/papers/www11.pdf}
     */
    constructor(elements, parameters = {}) {
        super(
            elements,
            // The neighborhood size is called `samples`; the default was previously keyed as `K`,
            // which nothing reads, so constructing without one left it `undefined`.
            /** @type {ParametersNNDescent} */ (
                Object.assign({ metric: euclidean, samples: 10, rho: 1, delta: 1e-3, seed: 1212 }, parameters)
            ),
        );
        this._N = elements.length;
        this._randomizer = new Randomizer(this._parameters.seed);
        this._sample_size = this._parameters.samples * this._parameters.rho;

        this._nndescent_elements = elements.map((e, i) => {
            return {
                value: e,
                index: i,
                flag: true,
            };
        });

        if (elements) {
            this.add(elements);
        }
    }

    /**
     * Samples Array A with sample size.
     *
     * @private
     * @template U
     * @param {U[]} A
     * @returns {U[]}
     */
    _sample(A) {
        const n = A.length;
        const sample_size = this._sample_size;
        if (sample_size > n) {
            return A;
        } else {
            const randomizer = this._randomizer;
            return randomizer.choice(A, sample_size);
        }
    }

    /**
     * Offers each of `u1`, `u2` to the other's neighborhood, under the distance between them.
     *
     * The distance has to be computed here. A neighbor object carries the distance to the point
     * whose list it came from, which says nothing about how far apart `u1` and `u2` are; reusing it
     * meant no new distance was ever evaluated and the heaps were ordered by unrelated values.
     *
     * @private
     * @param {KNNHeap<T>[]} B
     * @param {NNDescentNeighbor<T>} u1
     * @param {NNDescentNeighbor<T>} u2
     * @returns {number} How many neighborhoods changed.
     */
    _join(B, u1, u2) {
        if (u1.index === u2.index) return 0;
        const distance = this._parameters.metric(u1.value, u2.value);
        return (
            this._update(B[u1.index], { index: u2.index, value: u2.value, distance, flag: true }) +
            this._update(B[u2.index], { index: u1.index, value: u1.value, distance, flag: true })
        );
    }

    /**
     * @private
     * @param {KNNHeap<T>} B
     * @param {NNDescentNeighbor<T>} u
     * @returns {number}
     */
    _update(B, u) {
        if (B.set.has(u.index)) return 0;

        const worst = B.first;
        if (worst && B.length >= this._parameters.samples) {
            const dist = B._accessor(u);
            const worst_dist = B._accessor(worst.element);
            if (dist >= worst_dist) {
                return 0; // u is worse than the worst neighbor
            }
        }

        B.push(u);
        u.flag = true;
        if (B.length > this._parameters.samples) {
            B.pop();
        }
        return 1;
    }

    /**
     * @private
     * @param {(KNNHeap<T> | null)[]} B
     * @returns {NNDescentNeighbor<T>[][]}
     */
    _reverse(B) {
        const N = this._N;
        const R = new Array(N);
        for (let i = 0; i < N; i++) {
            R[i] = [];
        }
        for (let j = 0; j < N; j++) {
            const Bi = B[j];
            if (Bi) {
                const Bjdata = Bi.data();
                for (const neighbor of Bjdata) {
                    const v = neighbor.index;
                    R[v].push(neighbor);
                }
            }
        }
        return R;
    }

    /**
     * @param {T[]} elements
     * @returns {this}
     */
    add(elements) {
        const randomizer = this._randomizer;
        const metric = this._parameters.metric;
        const K = this._parameters.samples;
        const delta = this._parameters.delta;
        const N = elements.length;
        this._N = N;
        /** @type {KNNHeap<T>[]} */
        const B = [];
        this._B = B;
        // Hoisted: rebuilding this per point allocated N arrays of N objects.
        const indexed = elements.map((el, idx) => ({ el, idx }));
        for (let i = 0; i < N; i++) {
            const e = elements[i];
            // `flag: true` marks the entry as new. Without it nothing is ever eligible for the
            // local join below, so the graph stays at its random initialization.
            const sample = randomizer.choice(indexed, K).map((d) => {
                return { index: d.idx, distance: metric(d.el, e), value: d.el, flag: true };
            });
            const Bi = new KNNHeap(sample, (d) => d.distance, "max");
            B.push(Bi);
        }

        let c = Infinity;
        let old_c = -Infinity;
        while (c > delta * N * K && c !== old_c) {
            const old_ = new Array(N);
            const new_ = new Array(N);
            for (let i = 0; i < N; i++) {
                const Bi = B[i].data();
                const falseBs = Bi.filter((d) => !d.flag);
                const trueBs = this._sample(Bi.filter((d) => d.flag));
                for (const d of trueBs) {
                    d.flag = false;
                }
                old_[i] = new KNNHeap(falseBs, (d) => d.distance, "max");
                new_[i] = new KNNHeap(trueBs, (d) => d.distance, "max");
            }
            const old_reverse = this._reverse(old_);
            const new_reverse = this._reverse(new_);
            old_c = c;
            c = 0;
            for (let i = 0; i < N; i++) {
                for (const o of this._sample(old_reverse[i])) {
                    old_[i].push(o);
                }
                for (const n of this._sample(new_reverse[i])) {
                    new_[i].push(n);
                }

                const new_i = new_[i].data();
                const old_i = old_[i].data();
                const n1 = new_i.length;
                const n2 = old_i.length;
                for (let j = 0; j < n1; j++) {
                    const u1 = new_i[j];
                    // new x new is unordered, so start at j + 1 rather than joining each pair twice
                    for (let k = j + 1; k < n1; k++) {
                        c += this._join(B, u1, new_i[k]);
                    }
                    for (let k = 0; k < n2; k++) {
                        c += this._join(B, u1, old_i[k]);
                    }
                }
            }
        }
        this.nn = this._B.map((heap) => heap.data());
        return this;
    }

    /**
     * @param {T} x
     * @param {number} [k=5] Default is `5`
     * @returns {{ element: T, index: number; distance: number }[]}
     */
    search(x, k = 5) {
        const metric = this._parameters.metric;
        const N = this._N;
        const elements = this._elements;

        if (N === 0) return [];
        const xLength = x.length;

        // Initialize candidate pool
        const visited = new Set();
        /** @type {{index: number, dist: number, evaluated: boolean}[]} */
        let pool = [];

        // Randomly pick initial candidates
        const randomizer = this._randomizer;
        for (let i = 0; i < Math.min(N, Math.max(k, 10)); i++) {
            let rnd;
            do {
                rnd = randomizer.random_int % N;
            } while (visited.has(rnd));
            visited.add(rnd);

            const element = elements[rnd];
            if (!element || element.length !== xLength) continue;

            pool.push({
                index: rnd,
                dist: metric(x, element),
                evaluated: false,
            });
        }

        let searching = true;
        while (searching) {
            pool.sort((a, b) => a.dist - b.dist);
            // keep the top subset for exploration
            pool = pool.slice(0, Math.max(k * 2, 32));

            searching = false;
            // Bounded by the pool as it stands now: expansion appends to `pool`, and re-reading
            // `pool.length` here would walk those additions in the same pass, so the cap above
            // never applied and the walk degenerated into a scan of the whole graph.
            const frontier = pool.length;
            for (let i = 0; i < frontier; i++) {
                const candidate = pool[i];
                if (candidate.evaluated) continue;

                candidate.evaluated = true;
                searching = true;

                // get neighbors of this candidate from graph
                const neighbors = this.nn[candidate.index];
                if (!neighbors) continue;

                for (const neighbor of neighbors) {
                    const n_idx = neighbor.index;
                    if (!visited.has(n_idx)) {
                        visited.add(n_idx);
                        const element = elements[n_idx];
                        if (element && element.length === xLength) {
                            pool.push({
                                index: n_idx,
                                dist: metric(x, element),
                                evaluated: false,
                            });
                        }
                    }
                }
                // Don't break here! Look at more candidates per iteration for better convergence
                // break;
            }
        }

        pool.sort((a, b) => a.dist - b.dist);

        /** @type {{ element: T, index: number; distance: number }[]} */
        const result = [];
        for (let i = 0; i < Math.min(k, pool.length); i++) {
            const item = pool[i];
            result.push({
                element: elements[item.index],
                index: item.index,
                distance: item.dist,
            });
        }
        return result;
    }
}

/**
 * @template {number[] | Float64Array} U
 * @typedef {Object} HeapEntry
 * @property {NNDescentNeighbor<U>} element
 * @property {number} value
 */

/**
 * @template {number[] | Float64Array} U
 * @extends {Heap<NNDescentNeighbor<U>>}
 */
class KNNHeap extends Heap {
    /** @type {Set<number>} */
    set;

    /**
     * @param {NNDescentNeighbor<U>[]} elements
     * @param {(d: NNDescentNeighbor<U>) => number} accessor
     * @param {"max" | "min"} comparator
     */
    constructor(elements, accessor, comparator) {
        super(null, accessor, comparator);
        this.set = new Set();
        if (elements) {
            for (const element of elements) {
                this.push(element);
            }
        }
    }

    /**
     * @param {NNDescentNeighbor<U>} element
     * @returns {KNNHeap<U>}
     */
    push(element) {
        const set = this.set;
        if (set.has(element.index)) {
            return this;
        } else {
            set.add(element.index);
            super.push(element);
            return this;
        }
    }

    /** @returns {{ element: NNDescentNeighbor<U>; value: number } | null} */
    pop() {
        const result = super.pop();
        if (result?.element) {
            this.set.delete(result.element.index);
            return result;
        }
        return null;
    }

    /** @returns {NNDescentNeighbor<U>[]} */
    data() {
        return this._container.map((d) => d.element);
    }
}
