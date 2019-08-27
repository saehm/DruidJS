import { simultaneous_poweriteration } from "../linear_algebra/index";
import { norm } from "../matrix/index";
import { neumair_sum } from "../numerical/index";

/**
 * @class
 * @alias Matrix
 * @requires module:numerical/neumair_sum
 */
export class Matrix{
    /**
     * creates a new Matrix. Entries are stored in a Float64Array. 
     * @constructor
     * @memberof module:matrix
     * @alias Matrix
     * @param {number} rows - The amount of rows of the matrix.
     * @param {number} cols - The amount of columns of the matrix.
     * @param {(function|string|number)} value=0 - Can be a function with row and col as parameters, a number, or "zeros", "identity" or "I", or "center".
     *  - **function**: for each entry the function gets called with the parameters for the actual row and column.
     *  - **string**: allowed are
     *      - "zero", creates a zero matrix.
     *      - "identity" or "I", creates an identity matrix.
     *      - "center", creates an center matrix.
     *  - **number**: create a matrix filled with the given value.
     * @example
     * 
     * let A = new Matrix(10, 10, () => Math.random()); //creates a 10 times 10 random matrix.
     * let B = new Matrix(3, 3, "I"); // creates a 3 times 3 identity matrix.
     * @returns {Matrix} returns a {@link rows} times {@link cols} Matrix filled with {@link value}.
     */
    constructor(rows=null, cols=null, value=null) {
        this._rows = rows;
        this._cols = cols;
        this._data = null;
        if (rows && cols) {
            if (!value) {
                this._data = new Float64Array(rows * cols);
                return this;
            }
            if (typeof(value) === "function") {
                this._data = new Float64Array(rows * cols);
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        this._data[row * cols + col] = value(row, col);
                    }
                }
                return this;
            }
            if (typeof(value) === "string") {
                if (value === "zeros") {
                    return new Matrix(rows, cols, 0); 
                }
                if (value === "identity" || value === "I") {
                    this._data = new Float64Array(rows * cols);
                    for (let row = 0; row < rows; ++row) {
                        this._data[row * cols + row] = 1;
                    }
                    return this;
                }
                if (value === "center" && rows == cols) {
                    this._data = new Float64Array(rows * cols);
                    value = (i, j) => (i === j ? 1 : 0) - (1 / rows);
                    for (let row = 0; row < rows; ++row) {
                        for (let col = 0; col < cols; ++col) {
                            this._data[row * cols + col] = value(row, col);
                        }
                    }
                    return this;
                }
            }
            if (typeof(value) === "number") {
                this._data = new Float64Array(rows * cols);
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        this._data[row * cols + col] = value;
                    }
                }
                return this;
            }
        }
        return this;
        
    }

    /**
     * Creates a Matrix out of {@link A}.
     * @param {(Matrix|Array|Float64Array|number)} A - The matrix, array, or number, which should converted to a Matrix.
     * @param {string} [type = "row"] - If {@link A} is a Array or Float64Array, then type defines if it is a row- or a column vector. 
     * @returns {Matrix}
     * 
     * @example
     * let A = Matrix.from([[1, 0], [0, 1]]); //creates a two by two identity matrix.
     */
    static from(A, type="row") {
        if (A instanceof Matrix) {
            return A.clone();
        } else if (Array.isArray(A) || A instanceof Float64Array) {
            let m = A.length
            if (m === 0) throw "Array is empty";
            // 1d
            if (!Array.isArray(A[0]) && !(A[0] instanceof Float64Array)) {
                if (type === "row") {  
                    return new Matrix(1, m, (_, j) => A[j]);
                } else if (type === "col") {
                    return new Matrix(m, 1, (i) => A[i]);
                } else {
                    throw "1d array has NaN entries"
                }
            // 2d
            } else if (Array.isArray(A[0]) || A[0] instanceof Float64Array) {
                let n = A[0].length;
                for (let row = 0; row < m; ++row) {
                    if (A[row].length !== n) throw "various array lengths";
                }
                return new Matrix(m, n, (i, j) => A[i][j])
            }
        } else if (typeof(A) === "number") {
            return new Matrix(1, 1, A);
        } else {
            throw "error"
        }
    }

    /**
     * Returns the {@link row}th row from the Matrix.
     * @param {int} row 
     * @returns {Array}
     */
    row(row) {
        let result_row = new Array(this._cols);
        for (let col = 0; col < this._cols; ++col) {
            result_row[col] = this._data[row * this._cols + col];
        }
        return result_row;
    }

    /**
     * Sets the entries of {@link row}th row from the Matrix to the entries from {@link values}.
     * @param {int} row 
     * @param {Array} values 
     * @returns {Matrix}
     */
    set_row(row, values) {
        let cols = this._cols;
        if (Array.isArray(values) && values.length === cols) {
            let offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                this._data[offset + col] = values[col];
            }
        } else if (values instanceof Matrix && values.shape[1] === cols && values.shape[0] === 1) {
            let offset = row * cols;
            for (let col = 0; col < cols; ++col) {
                this._data[offset + col] = values._data[col];
            }
        }
        return this;
    }

    /**
     * Returns the {@link col}th column from the Matrix.
     * @param {int} col 
     * @returns {Array}
     */
    col(col) {
        let result_col = new Array(this._rows);
        for (let row = 0; row < this._rows; ++row) {
            result_col[row] = this._data[row * this._cols + col];
        }
        return result_col;
    }

    /**
     * Returns the {@link col}th entry from the {@link row}th row of the Matrix.
     * @param {int} row 
     * @param {int} col 
     * @returns {float64}
     */
    entry(row, col) {
        return this._data[row * this._cols + col];
    }

    /**
     * Sets the {@link col}th entry from the {@link row}th row of the Matrix to the given {@link value}.
     * @param {int} row 
     * @param {int} col 
     * @param {float64} value
     * @returns {Matrix}
     */
    set_entry(row, col, value) {
        this._data[row * this._cols + col] = value;
        return this;
    }

    /**
     * Returns a new transposed Matrix.
     * @returns {Matrix}
     */
    transpose() {
        let B = new Matrix(this._cols, this._rows, (row, col) => this.entry(col, row))
        return B;
    }

    /**
     * Returns a new transposed Matrix. Short-form of {@function transpose}.
     * @returns {Matrix}
     */
    get T() {
        return this.transpose();
    }

    /**
     * Returns the inverse of the Matrix.
     * @returns {Matrix}
     */
    inverse() {
        const rows = this._rows;
        const cols = this._cols;
        let B = new Matrix(rows, 2 * cols, (i,j) => {
            if (j >= cols) {
                return (i === (j - cols)) ? 1 : 0;
            } else {
                return this.entry(i, j);
            }
        });
        let h = 0; 
        let k = 0;
        while (h < rows && k < cols) {
            var i_max = 0;
            let max_val = -Infinity;
            for (let i = h; i < rows; ++i) {
                let val = Math.abs(B.entry(i,k));
                if (max_val < val) {
                    i_max = i;
                    max_val = val;
                }
            }
            if (B.entry(i_max, k) == 0) {
                k++;
            } else {
                // swap rows
                for (let j = 0; j < 2 * cols; ++j) {
                    let h_val = B.entry(h, j);
                    let i_val = B.entry(i_max, j);
                    B.set_entry(h, j, h_val);
                    B.set_entry(i_max, j, i_val);
                }
                for (let i = h + 1; i < rows; ++i) {
                    let f = B.entry(i, k) / B.entry(h, k);
                    B.set_entry(i, k, 0)
                    for (let j = k + 1; j < 2 * cols; ++j) {
                        B.set_entry(i, j, B.entry(i, j) - B.entry(h, j) * f);
                    }
                }
                h++;
                k++;
            }
        }

        for (let row = 0; row < rows; ++row) {
            let f = B.entry(row, row);
            for (let col = row; col < 2 * cols; ++col) {
                B.set_entry(row, col, B.entry(row, col) / f)
            }
        }
        
        for (let row = rows - 1; row >= 0; --row) {
            let B_row_row = B.entry(row, row);
            for (let i = 0; i < row; i++) {
                let B_i_row = B.entry(i, row);
                let f = B_i_row / B_row_row;
                for (let j = i; j < 2 * cols; ++j) {
                    let B_i_j = B.entry(i,j);
                    let B_row_j = B.entry(row, j);
                    B_i_j = B_i_j - B_row_j * f;
                    B.set_entry(i, j, B_i_j)
                }
            }
        }

        return new Matrix(rows, cols, (i,j) => B.entry(i, j + cols));
    }

    /**
     * Returns the dot product. If {@link B} is an Array or Float64Array then an Array gets returned. If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */
    dot(B) {
        if (B instanceof Matrix) {
            let A = this;
            if (A.shape[1] !== B.shape[0]) return undefined;
            let I = A.shape[1];
            let C = new Matrix(A.shape[0], B.shape[1], (row, col) => {
                let A_i = A.row(row);
                let B_i = B.col(col);
                for (let i = 0; i < I; ++i) {
                    A_i[i] *= B_i[i];
                }
                return neumair_sum(A_i);
            });
            return C;
        } else if (Array.isArray(B) || (B instanceof Float64Array)) {
            let rows = this._rows;
            if (B.length !== rows) return undefined;
            let C = new Array(rows)
            for (let row = 0; row < rows; ++row) {
                C[row] = neumair_sum(this.row(row).map(e => e * B[row]));
            }
            return C;
        } else {
            return undefined;
        }
    }

    /**
     * Computes the outer product from {@link this} and {@link B}.
     * @param {Matrix} B 
     * @returns {Matrix}
     */
    outer(B) {
        let A = this;
        let l = A._data.length;
        let r = B._data.length;
        if (l != r) return undefined;
        let C = new Matrix()
        C.shape = [l, l, (i, j) => {
            if (i <= j) {
                return A._data[i] * B._data[j];
            } else {
                return C.entry(j, i);
            }
        }]
        return C;
    }

    /**
     * Transforms A to bidiagonal form. (With Householder transformations)
     * @param {Matrix} A 
     */
    static bidiagonal(A) {
        A = A.clone()
        let [ rows, cols ] = A.shape;
        for (let k = 0; k < cols; ++k) {
            let x = A.col(k).slice(k);
            let x_norm = norm(x);
            let u = A.col(k).slice(k);
            u[0] = x[0] + (x[0] < 0 ? -1 : 1) * x_norm;
            let u_norm = norm(u);
            u = u.map(u_i => u_i /= u_norm);
            u = Matrix.from(u, "col");
            let A_ = A.get_block(k, k);
            A_ = A_.sub(u.mult(2).outer(u).dot(A_));
            A.set_block(k, k, A_); /** @todo repair numerical unstability? */
            // repair numerical unstability?
            for (let r = k + 1; r < rows; ++r) {
                A.set_entry(r, k, 0);
            }
            if (k <= cols - 2) {
                let y = A.row(k).slice(k + 1);
                let y_norm = norm(y);
                let v = A.row(k).slice(k + 1);
                v[0] = y[0] + (y[0] < 0 ? -1 : 1) * y_norm;
                let v_norm = norm(v)
                v = v.map(v_i => v_i /= v_norm);
                v = Matrix.from(v, "col");
                let _A = A.get_block(k, k + 1);
                _A = _A.sub(_A.dot(v.mult(2).outer(v))) 
                A.set_block(k, k + 1, _A);
                // repair numerical unstability?
                for (let c = k + 2; c < cols; ++c) {
                    A.set_entry(k, c, 0);
                }
            }
        }
        return A;
    }

    /**
     * Writes the entries of B in A at an offset position given by {@link offset_row} and {@link offset_col}.
     * @param {int} offset_row 
     * @param {int} offset_col 
     * @param {Matrix} B 
     * @returns {Matrix}
     */
    set_block(offset_row, offset_col, B) {
        let [ rows, cols ] = B.shape;
        for (let row = 0; row < rows; ++row) {
            if (row > this._rows) continue;
            for (let col = 0; col < cols; ++col) {
                if (col > this._cols) continue;
                this.set_entry(row + offset_row, col + offset_col, B.entry(row, col));
            }
        }
        return this;
    }

    /**
     * Extracts the entries of A at an offset position given by {@link offset_row} and {@link offset_col}.
     * @param {int} offset_row 
     * @param {int} offset_col 
     * @returns {Matrix}
     */
    get_block(offset_rows, offset_cols) {
        let [ rows, cols ] = this.shape;
        return new Matrix(rows - offset_rows, cols - offset_cols, (i, j) => this.entry(i + offset_rows, j + offset_cols));
    }

    /**
     * Applies a function to each entry of the matrix.
     * @param {function} f function takes 2 parameters, the value of the actual entry and a value given by the function {@link v}. The result of {@link f} gets writen to the Matrix.
     * @param {function} v function takes 2 parameters for row and col, and returns a value witch should be applied to the colth entry of the rowth row of the matrix.
     */
    _apply_array(f, v) {
        const data = this._data;
        const [ rows, cols ] = this.shape;
        for (let row = 0; row < rows; ++row) {
            const o = row * cols;
            for (let col = 0; col < cols; ++col) {
                const i = o + col;
                const d = data[i];
                data[i] = f(d, v(row, col));
            }
        }
        return this; 
    }

    _apply_rowwise_array(values, f) {
        return this._apply_array(f, (i, j) => values[j]);
    }

    _apply_colwise_array(values, f) {
        const data = this._data;
        const [ rows, cols ] = this.shape;
        for (let row = 0; row < rows; ++row) {
            const o = row * cols;
            for (let col = 0; col < cols; ++col) {
                const i = o + col;
                const d = data[i];
                data[i] = f(d, values[row]);
            }
        }
        return this; 
    }

    /*_apply_pointwise_number(value, f) {

    }*/

    _apply(value, f) {
        let data = this._data;
        if (value instanceof Matrix) {
            let [ value_rows, value_cols ] = value.shape;
            let [ rows, cols ] = this.shape;
            if (value_rows === 1) {
                if (cols !== value_cols) return undefined;
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(0, col))
                    }
                }
            } else if (value_cols === 1) {
                if (rows !== value_rows) return undefined;
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(row, 0))
                    }
                }
            } else if (rows == value_rows && cols == value_cols) {
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value.entry(row, col))
                    }
                }
            } else return undefined;
        } else if (Array.isArray(value)) {
            let rows = this._rows;
            let cols = this._cols;
            if (value.length === rows) {
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value[row]);
                    }
                }
            } else if (value.length === cols) {
                for (let row = 0; row < rows; ++row) {
                    for (let col = 0; col < cols; ++col) {
                        data[row * cols + col] = f(data[row * cols + col], value[col]);
                    }
                }
            } else {
                return undefined;
            }
        } else {
            for (let i = 0, n = this._rows * this._cols; i < n; ++i) {
                data[i] = f(data[i], value);
            }
        }
        return this;
    }

    /**
     * Clones the Matrix.
     * @returns {Matrix}
     */
    clone() {
        let B = new Matrix()
        B._rows = this._rows;
        B._cols = this._cols;
        B._data = this._data.slice(0);
        return B;
    }

    mult(value) {
        return this.clone()._apply(value, (a,b) => a * b)
    }

    divide(value) {
        return this.clone()._apply(value, (a,b) => a / b)
    }

    add(value) {
        return this.clone()._apply(value, (a,b) => a + b)
    }

    sub(value) {
        return this.clone()._apply(value, (a,b) => a - b)
    }

    /**
     * Returns the number of rows and columns of the Matrix.
     * @returns {Array} An Array in the form [rows, columns].
     */
    get shape() {
        return [this._rows, this._cols];
    }

    /**
     * Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.
     * @param {Array} parameter - takes an Array in the form [rows, cols, value], where rows and cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row and col) which has to return a value for the colth entry of the rowth row.
     * @returns {Matrix}
     */
    set shape([ rows, cols, value = () => 0 ]) {
        this._rows = rows;
        this._cols = cols;
        this._data = new Float64Array(rows * cols);
        for (let row = 0; row < rows; ++row) {
            for (let col = 0; col < cols; ++col) {
                this._data[row * cols + col] = value(row, col);
            }
        }
        return this;
    }

    /**
     * Returns the Matrix as a two-dimensional Array.
     * @returns {Array}
     */
    get to2dArray() {
        const rows = this._rows;
        const cols = this._cols;
        let result = new Array(rows)
        for (let row = 0; row < rows; ++row) {
            let result_col = new Array(cols)
            for (let col = 0; col < cols; ++col) {
                result_col[col] = this.entry(row, col);
            }
            result[row] = result_col;
        }
        return result;
    }

    /**
     * Returns the diagonal of the Matrix.
     * @returns {Array}
     */
    get diag() {
        const rows = this._rows;
        const cols = this._cols;
        const min_row_col = Math.min(rows, cols);
        let result = new Array(min_row_col)
        for (let i = 0; i < min_row_col; ++i) {
            result[i] = this.entry(i,i);
        }
        return result;
    }

    /**
     * Returns the mean of all entries of the Matrix.
     * @returns {float64}
     */
    get mean() {
        const data = this._data;
        const n = this._rows * this._cols;
        let sum = 0;
        for (let i = 0; i < n; ++i) {
            sum += data[i];
        }
        return sum / n;
    }

    /**
     * Returns the mean of each row of the matrix.
     * @returns {Array}
     */
    get meanRows() {
        const data = this._data;
        const rows = this._rows;
        const cols = this._cols;
        let result = [];
        for (let row = 0; row < rows; ++row) {
            result[row] = 0;
            for (let col = 0; col < cols; ++col) {
                result[row] += data[row * cols + col];
            }
            result[row] /= cols;
        }
        return result;
    }

    /** Returns the mean of each column of the matrix.
     * @returns {Array}
     */
    get meanCols() {
        const data = this._data;
        const rows = this._rows;
        const cols = this._cols;
        let result = [];
        for (let col = 0; col < cols; ++col) {
            result[col] = 0;
            for (let row = 0; row < rows; ++row) {
                result[col] += data[row * cols + col];
            }
            result[col] /= rows;
        }
        return result;
    }

    /**
     * Solves the equation {@link A}x = {@link b}. Returns the result x.
     * @param {Matrix} A - Matrix
     * @param {Matrix} b - Matrix
     * @returns {Matrix}
     */
    static solve(A, b) {
        let rows = A.shape[0];
        let { L: L, U: U } = Matrix.LU(A);//lu(A);
        let x = b.clone();
        
        // forward
        for (let row = 0; row < rows; ++row) {
            for (let col = 0; col < row - 1; ++col) {
                x.set_entry(0, row, x.entry(0, row) - L.entry(row, col) * x.entry(1, col));
            }
            x.set_entry(0, row, x.entry(0, row) / L.entry(row, row));
        }
        
        // backward
        for (let row = rows - 1; row >= 0; --row) {
            for (let col = rows - 1; col > row; --col) {
                x.set_entry(0, row, x.entry(0, row) - U.entry(row, col) * x.entry(0, col));
            }
            x.set_entry(0, row, x.entry(0, row) / U.entry(row, row));
        }

        return x;
    }

    /**
     * {@link L}{@link U} decomposition of the Matrix {@link A}. Creates two matrices, so that the dot product LU equals A.
     * @param {Matrix} A 
     * @returns {{L: Matrix, U: Matrix}} result - Returns the left triangle matrix {@link L} and the upper triangle matrix {@link U}.
     */
    static LU(A) {
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

    /**
     * Computes the {@link k} components of the SVD decomposition of the matrix {@link M}
     * @param {Matrix} M 
     * @param {int} [k=2] 
     * @returns {{U: Matrix, Sigma: Matrix, V: Matrix}}
     */
    static SVD(M, k=2) {
        let MtM = M.transpose().dot(M);
        let MMt = M.dot(M.transpose());
        let { eigenvectors: V, eigenvalues: Sigma } = simultaneous_poweriteration(MtM, k);
        let { eigenvectors: U } = simultaneous_poweriteration(MMt, k);

        return { U: U, Sigma: Sigma.map(sigma => Math.sqrt(sigma)), V: V };
        /*const [ rows, cols ] = matrix.shape;
        // 1.
        let U = new Matrix();
        U.shape = [rows, cols, (i, j) => i === j ? 1 : 0];
        let Sigma = matrix.clone();
        let Vt = new Matrix();
        Vt.shape = [cols, cols, (i, j) => i === j ? 1 : 0];
        let I = new Matrix();
        I.shape = [cols, cols, (i, j) => i === j ? 1 : 0];

        function householder_vector(x) {
            let dot_1on = norm(x.col(0).slice(1));
            let v = x.clone();
            v.set_entry(0, 0, 1);
            let beta = 0;
            if (dot_1on > 1e-63) {
                let x_0 = x.entry(0, 0);
                let x_norm = Math.sqrt(x_0 ** 2 + dot_1on);
                let v_0;
                if (x_0 <= 0) {
                    v_0 = x_0 - x_norm;
                } else {
                    v_0 = -dot_1on / (x_0 + x_norm);
                }
                v.set_entry(0, 0, v_0)
                beta = 2 * v_0 ** 2 / (dot_1on + v_0 ** 2)
                v.divide(v_0)
            }
            return {"v": v, "beta": beta}
        }

        function householder_matrix(size, col, v, beta) {
            let vvt = v.dot(v.transpose())
            let V = new Matrix()
            V.shape = [size, size, (r, c) => {
                if (r < col || c < col) {
                    return r === c ? 1 : 0;
                }
                return (r === c ? 1 : 0) - (beta * vvt.entry(r - col, c - col));
            }]

            console.log(V, vvt)
            return V;
        }

        for (let col = 0; col < cols; ++col) {
            let x = Matrix.from(Sigma.col(col), "col")
            let { v: u, beta: beta } = householder_vector(x);
            let Q = householder_matrix(cols, col, u, beta);
            console.log("u", u, Q)
            U = U.dot(Q)
            Sigma = Q.dot(Sigma);

            if (col < cols - 2) {
                let x = Matrix.from(Sigma.row(col).slice(1), "col")
                let { v: v, beta: beta } = householder_vector(x);
                Q = householder_matrix(cols, col + 1, v, beta);
                console.log("v", v, Q)
                Vt = Q.dot(Vt)
                Sigma = Sigma.dot(Q);
            }

            /*let u = Sigma.col(col);
            let u_norm = norm(u);
            let sign_u = u[col] >= 0 ? 1 : -1;
            u[col] = sign_u * (Math.abs(u[col]) + u_norm);
            u_norm = norm(u);
            //u = u.map(u_i => u_i / u_norm);
            for (let c = 0; c < cols; ++c) {
                u[c] = u[c] / u_norm;
            }
            u = Matrix.from(u, "col")
            let H = I.sub(u.dot(u.transpose()).mult(2))
            console.log(u.to2dArray, H.to2dArray)
            Sigma = H.dot(Sigma);
            U = U.dot(H.transpose())

            if (col < cols - 1) {

            }*/

            /*let u = Sigma.col(col).slice(col);
            let u_norm = norm(u);
            let sign_u = u[0] >= 0 ? 1 : -1;
            u[0] = sign_u * (Math.abs(u[0]) + u_norm);//+= sign_u * u_norm;
            u_norm = norm(u)//Math.sqrt(2 * u_norm * (u_norm + u[0]))//norm(u);
            u = u.map(u_i => u_i / u_norm);
            u = Matrix.from(u, "col")
            let A_k_m = new Matrix()
            A_k_m.shape = [rows - col, cols - col, (i, j) => Sigma.entry(i + col, j + col)];
            let U_k = u.dot(u.transpose().dot(A_k_m)).mult(2)
            for (let r = 0; r < rows - col; ++r) {
                for (let c = 0; c < cols - col; ++c) {
                    let val = Sigma.entry(col + r, col + c) - U_k.entry(r, c);
                    val = Math.abs(val) < 1e-15 ? 0 : val;
                    Sigma.set_entry(col + r, col + c, val)
                }
            }
            console.log(Sigma)



            /*if (col <= cols - 2) {
                let col_1 = col + 1
                let v = Sigma.row(col).slice(col_1);
                let v_norm = norm(v);
                let sign_v = v[0] >= 0 ? 1 : -1;
                v[0] += sign_v * (Math.abs(v[0]) + v_norm);//sign_v * v_norm;
                v_norm = norm(v);
                v = v.map(v_i => v_i / v_norm);
                v = Matrix.from(v, "col")
                let A_k1_n = new Matrix()
                A_k1_n.shape = [rows - col, cols - col_1, (i, j) => Sigma.entry(col + i, col_1 + j)];
                let V_k = A_k1_n.dot(v.dot(v.transpose())).mult(2)
                for (let r = 0; r < rows - col; ++r) {
                    for (let c = 0; c < cols - col_1; ++c) {
                        let val = Sigma.entry(col + r, col_1 + c) - V_k.entry(r, c)
                        console.log("V", [col + r, col_1 + c], 2 * Sigma.entry(col + r, col_1 + c) == V_k.entry(r, c))
                        val = Math.abs(val) < 1e-15 ? 0 : val;
                        Sigma.set_entry(col + r, col_1 + c, val)
                    }
                }
            }

            */
            /*console.log(Sigma)
        }*/
    }

    
}