import { euclidean } from "../metrics/index.js";
import { Randomizer } from "../util/index.js";
import { Heap } from "../datastructure/index.js";
/**
 * @class
 * @alias NNDescent
 */
export class NNDescent{
    /**
     * @constructor
     * @memberof module:knn
     * @alias NNDescent
     * @param {Array<*>=} elements - called V in paper.
     * @param {Function} [metric = euclidean] - called sigma in paper.
     * @param {Number} [K = 10] - number of neighbors {@link search} should return.
     * @param {Number} [rho = .8] - sample rate.
     * @param {Number} [delta = 0.0001] - precision parameter.
     * @param {Number} [seed = 1987] - seed for the random number generator.
     * @returns {NNDescent}
     * @see {@link http://www.cs.princeton.edu/cass/papers/www11.pdf}
     */
    constructor(elements, metric=euclidean, K = 10, rho = 1, delta = 1e-3, seed=19870307) {
        this._metric = metric;
        this._randomizer = new Randomizer(seed);
        this._K = K;
        this._rho = rho;
        this._sample_size = K * rho;
        this._delta = delta;
        if (elements) {
            this.add(elements)
        }
        return this;   
    }

    /**
     * Samples Array A with sample size.
     * @private
     * @param {Array<*>} A 
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
     * Updates the KNN heap and returns 1 if changed, or 0 if not.
     * @private
     * @param {KNNHeap} B 
     * @param {*} u 
     */
    _update(B, u) {
        if (B.push(u)) {
            u.flag = true;
            B.pop();
            return 1;
        } else {
            return 0;
        }
    }

    /**
     * Collects for each element where it is neighbor from.
     * @private
     * @param {Array<KNNHeap>} B 
     */
    _reverse(B) {
        const N = this._N;
        const R = new Array(N).fill().map(() => new Array());
        for (let i = 0; i < N; ++i) {
            for (let j = 0; j < N; ++j) {
                const Bjdata = B[j].data()
                const val = Bjdata.find(d => d.index === i)
                if (val) R[j].push(val);
            }
        }
        return R;
    }

    /**
     * 
     * @param {Array} elements 
     */
    add(elements) {
        this._elements = elements = elements.map((e, i) => {
            return {
                "element": e,
                "index": i,
                "flag": true,
            }
        })
        const randomizer = this._randomizer;
        const metric = this._metric;
        const K = this._K;
        const delta = this._delta;
        const N = this._N = elements.length;
        const B = this._B = new Array();
        // B[v] <-- Sample(V,K)
        for (let i = 0; i < N; ++i) {
            const e = elements[i];
            const sample = randomizer.choice(elements, K);
            const Bi = new KNNHeap(sample, (d) => metric(d.element, e.element), "max"); // "max" to pop the futherst elements away
            B.push(Bi);
        }

        // loop
        let c = Infinity;
        let old_c = -Infinity;
        //let min_iter = 10;
        //let max_iter = 20;
        //while (min_iter-- > 0 || (c < delta * N * K) && max_iter-- > 0) {
        while (c > (delta * N * K) && c != old_c) {
            // parallel for v e V do
            const old_ = new Array(N);
            const new_ = new Array(N);
            for (let i = 0; i < N; ++i) {
                const e = elements[i];
                const Bi = B[i].data();
                const falseBs = Bi.filter(d => !d.flag);
                const trueBs = this._sample(Bi.filter(d => d.flag));
                trueBs.forEach(d => d.flag = false);
                old_[i] = new KNNHeap(falseBs, (d) => metric(d.element, e.element), "max");
                new_[i] = new KNNHeap(trueBs, (d) => metric(d.element, e.element), "max");
            }
            const old_reverse = this._reverse(old_);
            const new_reverse = this._reverse(new_);
            old_c = c;
            c = 0;
            // parallel for v e V do
            for (let i = 0; i < N; ++i) {
                this._sample(old_reverse[i]).forEach(o => old_[i].push(o));
                this._sample(new_reverse[i]).forEach(n => new_[i].push(n));
                const new_i = new_[i].data();
                const old_i = old_[i].data();
                const n1 = new_i.length;
                const n2 = old_i.length;
                for (let j = 0; j < n1; ++j) {
                    const u1 = new_i[j];
                    const Bu1 = B[u1.index];
                    for (let k = 0; k < n1; ++k) {
                        const u2 = new_i[k];
                        if (u1 == u2) continue;
                        const Bu2 = B[u2.index];
                        c += this._update(Bu2, u1);
                        c += this._update(Bu1, u2);
                    }
                    for (let k = 0; k < n2; ++k) {
                        const u2 = old_i[k];
                        if (u1 == u2) continue;
                        const Bu2 = B[u2.index];
                        c += this._update(Bu2, u1);
                        c += this._update(Bu1, u2);
                    }
                }
            }
        } 
        return this;
    }

    /**
     * @todo not implemented yet
     * @param {*} x 
     * @param {*} k 
     */
    search(x, k=5) {
        return this._B[this._randomizer.random_int % (this._N - 1)].toArray().slice(0, k);
    }

    search_index(i, k=5) {
        const B = this._B[i];
        const result = B.raw_data().sort((a, b) => a.value - b.value).slice(-k);
        return result;
    }
}

class KNNHeap extends Heap{
    constructor(elements, accessor, comparator) {
        super(null, accessor, comparator)
        this.set = new Set();
        if (elements) {
            for (const element of elements) {
                this.push(element);
            }
        }
    }

    push(element) {
        const set = this.set;
        if (set.has(element)){
            return false;
        } else {
            set.add(element);
            super.push(element);
            return true;
        }    
    }

    pop() {
        super.pop().element;
        //const element = super.pop().element;
        // once popped it should not return into the heap.
        // used as max heap. therefore, if popped the furthest 
        // element in the knn list gets removed.
        // this.set.delete(element); 
    }
}