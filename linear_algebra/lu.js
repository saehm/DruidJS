import { Matrix, norm, linspace } from "../matrix/index.js";
import { euclidean } from "../metrics/index.js";
import { neumair_sum } from "../numerical/index.js";

// crout algorithm
// https://en.wikipedia.org/wiki/Crout_matrix_decomposition
export default function(A) {
    let rows = A.shape[0];
    let L = new Matrix(rows, rows, "zeros");
    let U = new Matrix(rows, rows, "identity");
    let sum;

    for (let j = 0; j < rows; ++j) {
        for (let i = j; i < rows; ++i) {
            sum = 0
            for (let k = 0; k < j; ++k) {
                sum += L.entry(i, k) * U.entry(k, j)
            }
            /*sum = neumair_sum(linspace(0, j).map((k) => L.entry(i, k) * U.entry(k, j)))
            console.log("lsum1", sum)
            sum = []
            for (let k = 0; k < j; ++k) {
                sum.push(L.entry(i, k) * U.entry(k, j))
            }
            sum = neumair_sum(sum)
            console.log("lsum2", sum)*/
            L.set_entry(i, j, A.entry(i, j) - sum);
        }
        for (let i = j; i < rows; ++i) {
            if (L.entry(j, j) === 0) {
                return undefined;
            }
            sum = 0
            for (let k = 0; k < j; ++k) {
                sum += L.entry(j, k) * U.entry(k, i)
            }
            /*sum = neumair_sum(linspace(0, j).map((k) => L.entry(j, k) * U.entry(k, i)))
            console.log("usum1", sum)
            sum = []
            for (let k = 0; k < j; ++k) {
                sum.push(L.entry(j, k) * U.entry(k, i))
            }
            sum = neumair_sum("usum2", sum)
            console.log(sum)*/
            U.set_entry(j, i, (A.entry(j, i) - sum) / L.entry(j, j));
        }
    }

    return { L: L, U: U };
}

// doolittle algorithm
/*export default function(M) {
    let [rows, cols] = M.shape;
    let P = new Matrix();
    P.shape = [rows + 1, 1, i => i];
    let A = M.clone();
    let L = new Matrix();
    let U = new Matrix();
    let I = new Matrix();
    I.shape = [rows, rows, (i, j) => i === j ? 1 : 0];

    for (let i = 0; i < rows; ++i) {
        let max_A = 0;
        let i_max = i;
        let abs_A;
        for (let k = i; k < rows; ++k) {
            abs_A = Math.abs(A.entry(k, i))
            if (abs_A > max_A) {
                [ max_A, i_max ] = [ abs_A, k ];
            }

            if (max_A < 1e-12) return undefined;

            if (i_max !== i) {
                let p = P.row(i);
                P.set_row(i, P.row(i_max));
                P.set_row(i_max, p);

                let A_row = A.row(i);
                A.set_row(i, A.row(i_max));
                A.set_row(i_max, A_row);

                P.set_entry(rows + 1, 0, P.entry(rows + 1, 0) + 1)
            }
        }

        for (let j = i + 1; j < rows; ++j) {
            A.set_entry(j, i,  A.entry(j, i) / A.entry(i, i));
            for (let k = i + 1; k < rows; ++k) {
                A.set_entry(j, k, A.entry(j, k) - A.entry(j, i) * A.entry(i, k));
            }
        }
    }

    L.shape = [rows, rows, (i, j) => {
        if ( i > j ) return A.entry(i, j);
        if ( i === j ) return A.entry(i, j) + 1;
        return 0
    }]

    U = A.sub(L)

    return {L: L, U: U, P: P};   
}*/

