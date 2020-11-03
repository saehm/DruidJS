import { distance_matrix as dmatrix} from "../matrix/index";
import { euclidean } from '../metrics/index';

/**
 * 
 * @param {*} A 
 * @param {*} k 
 * @param {*} distance_matrix 
 * @param {*} metric 
 */
export default function(A, k, distance_matrix = null, metric = euclidean) {
    const rows = A.shape[0];
    let D = distance_matrix ?? dmatrix(A, metric)
    /* for (let i = 0; i < n; ++i) {
        D[i] = Array.from(D[i]).map((_,j) => {
                return {
                    i: i, j: j, distance: D[i][j]
                }
            })
            .sort((a, b) => a.distance - b.distance)
            .slice(1, k + 1)
    } */
    let nN = new Array(rows);
    for (let row = 0; row < rows; ++row) {
        nN[row] = Array.from(D.row(row)).map((distance, col) => {
                return {
                    "i": row,
                    "j": col,
                    "distance": distance,
                }
            })
            .sort((a, b) => a.distance - b.distance)
            .slice(1, k + 1);
    }
    return nN;
} 
