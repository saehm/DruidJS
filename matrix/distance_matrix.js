import { euclidean } from "../metrics/index"

export default function(A, metric = euclidean) {
    if (metric === undefined) return undefined;
    let n = A.length;
    let D = new Array(n);
    for (let i = 0; i < n; ++i) {
        D[i] = new Array(n);
    }
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            D[i][j] = D[j][i] = metric(A[i], A[j]);
        }
    }
    return D;
}
