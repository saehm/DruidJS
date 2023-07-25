import { Heap } from "../datastructure/index.js";

/**
 * 
 * @param {Array} A array of objects [{index: other_point, distance: distance to other_point} ** k]
 */
export default function(A) {
    console.log(A)
    let result = [];
    /* for (let i = 0, n = A.length; i < n; ++i) {
        let Q = new Heap(null, d => d.distance, "min");
        let source = A[i];
        let dist = {};
        let prev = {};
        for (let j = 0, m = source.length; j < m; ++j) {
            let v = source[j].index;
            dist[v] = Infinity;
            prev[v] = undefined;
            Q.push(source[j])
        }

        while (!Q.empty) {
            let u = Q.pop();

        }

        
    } */
    
}