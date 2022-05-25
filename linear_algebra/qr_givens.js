import { Matrix } from "../matrix/index.js";

/**
 * Computes the QR Decomposition of the Matrix {@link A} with givens rotation.
 * @memberof module:linear_algebra
 * @alias qr_givens
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_Givens_rotations}
 */
export default function (A) {
    const [rows, cols] = A.shape;
    let Q = new Matrix(rows, rows, "I");
    let R = A.clone();

    for (let j = 0; j < cols; ++j) {
        for (let i = rows - 1; i > j; --i) {
            const [c, s] = givensrotation(R.entry(i - 1, j), R.entry(i, j));
            if (c == 1 && s == 0) continue;
            const Gji = new Matrix(rows, rows, "I");
            Gji.set_entry(i - 1, i - 1, c);
            Gji.set_entry(i - 1, i, -s);
            Gji.set_entry(i, i - 1, s);
            Gji.set_entry(i, i, c);
            R = Gji.transDot(R);
            Q = Q.dot(Gji);
        }
    }
    return { R, Q };
}

function givensrotation(a, b) {
    if (b == 0) {
        return [1, 0];
    } else {
        if (Math.abs(b) > Math.abs(a)) {
            const r = a / b;
            const s = 1 / Math.sqrt(1 + r ** 2);
            const c = s * r;
            return [c, s];
        } else {
            const r = b / a;
            const c = 1 / Math.sqrt(1 + r ** 2);
            const s = c * r;
            return [c, s];
        }
    }
}
