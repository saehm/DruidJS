import { euclidean } from "../metrics/index";
import { Randomizer } from "../util/index";
import { Heap } from "../datastructure/index";
/**
 * @memberof module:knn
 */
export class NNDescent{
    /**
     * @see {@link http://www.cs.princeton.edu/cass/papers/www11.pdf}
     * @see {@link https://github.com/lmcinnes/pynndescent/blob/master/pynndescent/pynndescent_.py}
     */
    constructor(elements=null, metric=euclidean, seed=19870307) {
        this._metric = metric;
        const rand = this._randomizer = new Randomizer(seed)
        if (elements) {
            const k = this._k = 5;
            const B = this._knn_list = new Heap(rand.choice(elements, k), d => d, "min");
            let c = true;
            while (c) {
                let R = this._reverse(B);
                c = false
            }

        }
        

        return this;   
    }

    _reverse(B) {
        return B;
    }

    /**
     * 
     * @param {Array} elements 
     */
    add(elements) {
        elements;
        return this;
    }

    search(x, k=5) {
        return new Heap(x, k)
    }
}