import { euclidean } from "../metrics/index"

export default function(A, metric = euclidean) {
    let distance = metric;
    if (distance === undefined) return undefined;
    let n = A.length;
    let D = new Array(n);
    for (let i = 0; i < n; ++i) {
        D[i] = new Array(n);
    }
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            D[i][j] = D[j][i] = distance(A[i], A[j]);
        }
    }
    return D;
}
