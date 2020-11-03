import { euclidean } from "../metrics/index"
import { Matrix } from "./Matrix";

export default function(A, metric = euclidean) {
    let n = A.shape[0];
    /* let D = new Array(n);
    for (let i = 0; i < n; ++i) {
        D[i] = new Float64Array(n);
    }
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            D[i][j] = D[j][i] = metric(A[i], A[j]);
        }
    } */
    const D = new Matrix(n, n);
    for (let i = 0; i < n; ++i) {
        const A_i = A.row(i);
        for (let j = i + 1; j < n; ++j) {
            const dist = metric(A_i, A.row(j));
            D.set_entry(i, j, dist);
            D.set_entry(j, i, dist);
        }
    }
    return D;
}
