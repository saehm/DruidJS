import { euclidean } from "../metrics/index.js";
import { Matrix, norm } from "../matrix/index.js";

export default function(A, k = 2, max_iterations = 100, metric = euclidean) {
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    let n = A.shape[0]
    let v = new Matrix(n, 1, () => Math.random());
    v = v.divide(norm(v._data, metric))
    
    let w_ = A.dot(v);
    let a = w_.dot(v);
    let w = w_.sub(a.dot(v));

    while (max_iterations--) {
        let b = norm(w._data, euclidean);
        let v;
        if (b === 0) {
            let v_rand
        } 
    }
    return {
        "eigenvector" : r, 
        "eigenvalue": u / l
    };
}

