import { qr } from "./index";
import { Matrix } from "../matrix/Matrix";
import { Randomizer } from "../util/index";
import { neumair_sum } from "../numerical/index";

export default function(A, k = 2, max_iterations = 100, seed = 19870307) {
    let randomizer = new Randomizer(seed)
    if (!(A instanceof Matrix)) A = Matrix.from(A);
    let n = A.shape[0]
    let { Q: Q, R: R } = qr(new Matrix(n, k, () => randomizer.random));
    while(max_iterations--) {
        let oldR = R.clone();
        let Z = A.dot(Q);
        let QR = qr(Z);
        [ Q, R ] = [ QR.Q, QR.R ]; 
        if (neumair_sum(R.sub(oldR).diag) / n < 1e-12) {
            max_iterations = 0;
        }        
    }

    let eigenvalues = R.diag;
    let eigenvectors = Q.transpose().to2dArray;//.map((d,i) => d.map(dd => dd * eigenvalues[i]))
    return {
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors
    };
}

