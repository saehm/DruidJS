import { distance_matrix as dmatrix} from "../matrix/index";
import { euclidean } from '../metrics/index';

export default function(A, k, distance_matrix = null, metric = euclidean) {
    let n = A.length
    let D = distance_matrix || dmatrix(A, metric)
    for (let i = 0; i < n; ++i) {
        D[i] = D[i].map((d,j) => {
            return {
                i: i, j: j, distance: D[i][j]
            }
        }).sort((a, b) => a.distance - b.distance)
        .slice(1, k + 1)
    }
    return D
} 
