import { simultaneous_poweriteration} from "../linear_algebra/index";
import { Matrix } from "../matrix/index";

export class PCA{
    constructor(X, d=2) {
        this.X = X;
        this.d = d;
    }

    transform() {
        let X = this.X;
        let D = X.shape[1];
        let O = new Matrix(D, D, "center");
        let X_cent = X.dot(O);

        let C = X_cent.transpose().dot(X_cent)
        let { eigenvectors: V } = simultaneous_poweriteration(C, this.d)
        console.log(V)
        V = Matrix.from(V).transpose()
        this.Y = X.dot(V)
        return this.Y
    }

    get projection() {
        return this.Y
    }

    

} 