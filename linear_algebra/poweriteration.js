import { default as euclidean } from "../metrics/euclidean"
import { default as random_array } from "../matrix/random_array"
import { divide, dot, norm, transpose } from "../matrix/index"

export default function(A, max_iterations = 100, metric = euclidean) {
    let n = A.length;
    let r = random_array(n,0);
    r = divide(r, norm(r, metric))
    while (max_iterations--) {
        let r_next = dot(A, r);
        r = divide(r_next, norm(r_next, metric));
    }
    let u = dot(dot(transpose(A), r), r)
    let l = dot(r, r);
    return {
        "eigenvector" : r, 
        "eigenvalue": u / l
    };
}

