// https://renecutura.eu v0.3.16 Copyright 2021 Rene Cutura
/**
 * Computes the euclidean distance (l_2) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias euclidean
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the euclidean distance between {@link a} and {@link b}.  
 */
function euclidean(t,e){return Math.sqrt(euclidean_squared(t,e))}
/**
 * Numerical stable summation with the Kahan summation algorithm.
 * @memberof module:numerical
 * @alias kahan_sum
 * @param {Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm}
 */function kahan_sum(t){let e,r,s=t.length,i=0,n=0;for(let o=0;o<s;++o)e=t[o]-n,r=i+e,n=r-i-e,i=r;return i}
/**
 * Numerical stable summation with the Neumair summation algorithm.
 * @memberof module:numerical
 * @alias neumair_sum
 * @param {Array} summands - Array of values to sum up.
 * @returns {number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements}
 */function neumair_sum(t){let e=t.length,r=0,s=0;for(let i=0;i<e;++i){let e=t[i],n=r+e;Math.abs(r)>=Math.abs(e)?s+=r-n+e:s+=e-n+r,r=n}return r+s}
/**
 * Computes the squared euclidean distance (l_2^2) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias euclidean_squared
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the squared euclidean distance between {@link a} and {@link b}.  
 */function euclidean_squared(t,e){if(t.length!=e.length)return;let r=t.length,s=new Array(r);for(let i=0;i<r;++i){let r=t[i],n=e[i];s[i]=(r-n)*(r-n)}return neumair_sum(s)}
/**
 * Computes the cosine distance (not similarity) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias cosine
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @example
 * druid.cosine([1,0],[1,1]) == 0.7853981633974484 == π/4
 * @returns {Number} The cosine distance between {@link a} and {@link b}.
 */function cosine(t,e){if(t.length!==e.length)return;let r=t.length,s=0,i=0,n=0;for(let o=0;o<r;++o)s+=t[o]*e[o],i+=t[o]*t[o],n+=e[o]*e[o];return Math.acos(s/(Math.sqrt(i)*Math.sqrt(n)))}
/**
 * Computes the manhattan distance (l_1) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias manhattan
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the manhattan distance between {@link a} and {@link b}.  
 */function manhattan(t,e){if(t.length!=e.length)return;let r=t.length,s=0;for(let i=0;i<r;++i)s+=Math.abs(t[i]-e[i]);return s}
/**
 * Computes the chebyshev distance (l_∞) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias chebyshev
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the chebyshev distance between {@link a} and {@link b}.  
 */function chebyshev(t,e){if(t.length!=e.length)return;let r=t.length,s=[];for(let i=0;i<r;++i)s.push(Math.abs(t[i]-e[i]));return Math.max(...s)}
/**
 * Computes the canberra distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias canberra
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} The canberra distance between {@link a} and {@link b}.
 * @see {@link https://en.wikipedia.org/wiki/Canberra_distance}
 */function canberra(t,e){if(t.length!==e.length)return;let r=t.length,s=0;for(let i=0;i<r;++i)s+=Math.abs(t[i]-e[i])/(Math.abs(t[i])+Math.abs(e[i]));return s}
/**
 * Computes the jaccard distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias jaccard
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the jaccard distance between {@link a} and {@link b}.  
 */function jaccard(t,e){if(t.length!=e.length)return;const r=t.length;let s=0,i=0;for(let n=0;n<r;++n){const r=0!=t[n],o=0!=e[n];s+=r||o,i+=r&&o}return(s-i)/s}
/**
 * Computes the hamming distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias hamming
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the hamming distance between {@link a} and {@link b}.  
 */function hamming(t,e){if(t.length!=e.length)return;const r=t.length;let s=0;for(let i=0;i<r;++i){s+=t[i]!=e[i]}return s/r}
/**
 * Computes the Sokal-Michener distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias sokal_michener
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the Sokal-Michener distance between {@link a} and {@link b}.  
 */function sokal_michener(t,e){if(t.length!=e.length)return;const r=t.length;let s=0;for(let i=0;i<r;++i){s+=0!=t[i]!=(0!=e[i])}return 2*s/(r+s)}
/**
 * Computes the yule distance between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias yule
 * @param {Array<Number>} a 
 * @param {Array<Number>} b 
 * @returns {Number} the yule distance between {@link a} and {@link b}.  
 */function yule(t,e){if(t.length!=e.length)return;const r=t.length;let s=0,i=0,n=0;for(let o=0;o<r;++o){const r=0!=t[o],a=0!=e[o];s+=r&&a,i+=r&&!a,n+=!r&&r}return 0==i||0==n?0:2*i*n/(s*(r-s-i-n)+i*n)}
/**
 * 
 * @param {*} A 
 * @param {*} k 
 * @param {*} distance_matrix 
 * @param {*} metric 
 */function k_nearest_neighbors(t,e,r=null,s=euclidean){const i=t.shape[0];let n=r??distance_matrix(t,s),o=new Array(i);
/* for (let i = 0; i < n; ++i) {
        D[i] = Array.from(D[i]).map((_,j) => {
                return {
                    i: i, j: j, distance: D[i][j]
                }
            })
            .sort((a, b) => a.distance - b.distance)
            .slice(1, k + 1)
    } */for(let t=0;t<i;++t)o[t]=Array.from(n.row(t)).map(((e,r)=>({i:t,j:r,distance:e}))).sort(((t,e)=>t.distance-e.distance)).slice(1,e+1);return o}
/**
 * @class
 * @alias Matrix
 * @requires module:numerical/neumair_sum
 */class Matrix{
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
constructor(t=null,e=null,r=null){if(this._rows=t,this._cols=e,this._data=null,t&&e){if(!r)return this._data=new Float64Array(t*e),this;if("function"==typeof r){this._data=new Float64Array(t*e);for(let s=0;s<t;++s)for(let t=0;t<e;++t)this._data[s*e+t]=r(s,t);return this}if("string"==typeof r){if("zeros"===r)return new Matrix(t,e,0);if("identity"===r||"I"===r){this._data=new Float64Array(t*e);for(let r=0;r<t;++r)this._data[r*e+r]=1;return this}if("center"===r&&t==e){this._data=new Float64Array(t*e),r=(e,r)=>(e===r?1:0)-1/t;for(let s=0;s<t;++s)for(let t=0;t<e;++t)this._data[s*e+t]=r(s,t);return this}}if("number"==typeof r){this._data=new Float64Array(t*e);for(let s=0;s<t;++s)for(let t=0;t<e;++t)this._data[s*e+t]=r;return this}}return this}
/**
     * Creates a Matrix out of {@link A}.
     * @param {(Matrix|Array|Float64Array|number)} A - The matrix, array, or number, which should converted to a Matrix.
     * @param {"row"|"col"|"diag"} [type = "row"] - If {@link A} is a Array or Float64Array, then type defines if it is a row- or a column vector. 
     * @returns {Matrix}
     * 
     * @example
     * let A = Matrix.from([[1, 0], [0, 1]]); //creates a two by two identity matrix.
     * let S = Matrix.from([1, 2, 3], "diag"); // creates a three by three matrix with 1, 2, 3 on its diagonal.
     */static from(t,e="row"){if(t instanceof Matrix)return t.clone();if(!(Array.isArray(t)||t instanceof Float64Array)){if("number"==typeof t)return new Matrix(1,1,t);throw"error"}{let r=t.length;if(0===r)throw"Array is empty";
// 1d
if(!(Array.isArray(t[0])||t[0]instanceof Float64Array)){if("row"===e)return new Matrix(1,r,((e,r)=>t[r]));
// 2d
if("col"===e)return new Matrix(r,1,(e=>t[e]));if("diag"===e)return new Matrix(r,r,((e,r)=>e==r?t[e]:0));throw"1d array has NaN entries"}if(Array.isArray(t[0])||t[0]instanceof Float64Array){let e=t[0].length;for(let s=0;s<r;++s)if(t[s].length!==e)throw"various array lengths";return new Matrix(r,e,((e,r)=>t[e][r]))}}}
/**
     * Returns the {@link row}th row from the Matrix.
     * @param {int} row 
     * @returns {Array}
     */row(t){
/* let result_row = new Array(this._cols);
        for (let col = 0; col < this._cols; ++col) {
            result_row[col] = this._data[row * this._cols + col];
        }
        return result_row; */
const e=this._data,r=this._cols;return e.subarray(t*r,(t+1)*r)}
/**
     * Returns an generator yielding each row of the Matrix.
     */*iterate_rows(){const t=this._cols,e=this._rows,r=this._data;for(let s=0;s<e;++s)yield r.subarray(s*t,(s+1)*t)}
/**
     * Makes a {@link Matrix} object an iterable object.
     */*[Symbol.iterator](){for(const t of this.iterate_rows())yield t}
/**
     * Sets the entries of {@link row}th row from the Matrix to the entries from {@link values}.
     * @param {int} row 
     * @param {Array} values 
     * @returns {Matrix}
     */set_row(t,e){let r=this._cols;if(Array.isArray(e)&&e.length===r){let s=t*r;for(let t=0;t<r;++t)this._data[s+t]=e[t]}else if(e instanceof Matrix&&e.shape[1]===r&&1===e.shape[0]){let s=t*r;for(let t=0;t<r;++t)this._data[s+t]=e._data[t]}return this}
/**
     * Returns the {@link col}th column from the Matrix.
     * @param {int} col 
     * @returns {Array}
     */col(t){let e=new Float64Array(this._rows);for(let r=0;r<this._rows;++r)e[r]=this._data[r*this._cols+t];return e}
/**
     * Returns the {@link col}th entry from the {@link row}th row of the Matrix.
     * @param {int} row 
     * @param {int} col 
     * @returns {float64}
     */entry(t,e){return this._data[t*this._cols+e]}
/**
     * Sets the {@link col}th entry from the {@link row}th row of the Matrix to the given {@link value}.
     * @param {int} row 
     * @param {int} col 
     * @param {float64} value
     * @returns {Matrix}
     */set_entry(t,e,r){return this._data[t*this._cols+e]=r,this}
/**
     * Returns a new transposed Matrix.
     * @returns {Matrix}
     */transpose(){return new Matrix(this._cols,this._rows,((t,e)=>this.entry(e,t)))}
/**
     * Returns a new transposed Matrix. Short-form of {@function transpose}.
     * @returns {Matrix}
     */get T(){return this.transpose()}
/**
     * Returns the inverse of the Matrix.
     * @returns {Matrix}
     */inverse(){const t=this._rows,e=this._cols;let r=new Matrix(t,2*e,((t,r)=>r>=e?t===r-e?1:0:this.entry(t,r))),s=0,i=0;for(;s<t&&i<e;){var n=0;let o=-1/0;for(let e=s;e<t;++e){let t=Math.abs(r.entry(e,i));o<t&&(n=e,o=t)}if(0==r.entry(n,i))i++;else{
// swap rows
for(let t=0;t<2*e;++t){let e=r.entry(s,t),i=r.entry(n,t);r.set_entry(s,t,e),r.set_entry(n,t,i)}for(let n=s+1;n<t;++n){let t=r.entry(n,i)/r.entry(s,i);r.set_entry(n,i,0);for(let o=i+1;o<2*e;++o)r.set_entry(n,o,r.entry(n,o)-r.entry(s,o)*t)}s++,i++}}for(let s=0;s<t;++s){let t=r.entry(s,s);for(let i=s;i<2*e;++i)r.set_entry(s,i,r.entry(s,i)/t)}for(let s=t-1;s>=0;--s){let t=r.entry(s,s);for(let i=0;i<s;i++){let n=r.entry(i,s)/t;for(let t=i;t<2*e;++t){let e=r.entry(i,t);e-=r.entry(s,t)*n,r.set_entry(i,t,e)}}}return new Matrix(t,e,((t,s)=>r.entry(t,s+e)))}
/**
     * Returns the dot product. If {@link B} is an Array or Float64Array then an Array gets returned. If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */dot(t){if(t instanceof Matrix){let e=this;if(e.shape[1]!==t.shape[0])throw`A.dot(B): A is a ${e.shape.join(" x ")}-Matrix, B is a ${t.shape.join(" x ")}-Matrix: \n                A has ${e.shape[1]} cols and B ${t.shape[0]} rows. \n                Must be equal!`;let r=e.shape[1];return new Matrix(e.shape[0],t.shape[1],((s,i)=>{const n=e.row(s),o=t.col(i);let a=0;for(let t=0;t<r;++t)a+=n[t]*o[t];return a}))}if(Array.isArray(t)||t instanceof Float64Array){let e=this._rows;if(t.length!==e)throw`A.dot(B): A has ${e} cols and B has ${t.length} rows. Must be equal!`;let r=new Array(e);for(let s=0;s<e;++s)r[s]=neumair_sum(this.row(s).map((e=>e*t[s])));return r}throw"B must be Matrix or Array"}
/**
     * Computes the outer product from {@link this} and {@link B}.
     * @param {Matrix} B 
     * @returns {Matrix}
     */outer(t){let e=this,r=e._data.length;if(r!=t._data.length)return;let s=new Matrix;return s.shape=[r,r,(r,i)=>r<=i?e._data[r]*t._data[i]:s.entry(i,r)],s}
/**
     * Appends matrix {@link B} to the matrix.
     * @param {Matrix} B - matrix to append.
     * @param {"horizontal"|"vertical"|"diag"} [type = "horizontal"] - type of concatenation.
     * @returns {Matrix}
     * @example
     * 
     * let A = Matrix.from([[1, 1], [1, 1]]); // 2 by 2 matrix filled with ones.
     * let B = Matrix.from([[2, 2], [2, 2]]); // 2 by 2 matrix filled with twos.
     * 
     * A.concat(B, "horizontal"); // 2 by 4 matrix. [[1, 1, 2, 2], [1, 1, 2, 2]]
     * A.concat(B, "vertical"); // 4 by 2 matrix. [[1, 1], [1, 1], [2, 2], [2, 2]]
     * A.concat(B, "diag"); // 4 by 4 matrix. [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2], [0, 0, 2, 2]]
     */concat(t,e="horizontal"){const r=this,[s,i]=r.shape,[n,o]=t.shape;if("horizontal"==e){if(s!=n)throw`A.concat(B, "horizontal"): A and B need same number of rows, A has ${s} rows, B has ${n} rows.`;const e=new Matrix(s,i+o,"zeros");return e.set_block(0,0,r),e.set_block(0,i,t),e}if("vertical"==e){if(i!=o)throw`A.concat(B, "vertical"): A and B need same number of columns, A has ${i} columns, B has ${o} columns.`;const e=new Matrix(s+n,i,"zeros");return e.set_block(0,0,r),e.set_block(s,0,t),e}if("diag"==e){const e=new Matrix(s+n,i+o,"zeros");return e.set_block(0,0,r),e.set_block(s,i,t),e}throw`type must be "horizontal" or "vertical", but type is ${e}!`}
/**
     * Writes the entries of B in A at an offset position given by {@link offset_row} and {@link offset_col}.
     * @param {int} offset_row 
     * @param {int} offset_col 
     * @param {Matrix} B 
     * @returns {Matrix}
     */set_block(t,e,r){let[s,i]=r.shape;for(let n=0;n<s;++n)if(!(n>this._rows))for(let s=0;s<i;++s)s>this._cols||this.set_entry(n+t,s+e,r.entry(n,s));return this}
/**
     * Extracts the entries from the {@link start_row}th row to the {@link end_row}th row, the {@link start_col}th column to the {@link end_col}th column of the matrix.
     * If {@link end_row} or {@link end_col} is empty, the respective value is set to {@link this.rows} or {@link this.cols}.
     * @param {Number} start_row 
     * @param {Number} start_col
     * @param {Number} [end_row = null]
     * @param {Number} [end_col = null] 
     * @returns {Matrix} Returns a end_row - start_row times end_col - start_col matrix, with respective entries from the matrix.
     * @example
     * 
     * let A = Matrix.from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]); // a 3 by 3 matrix.
     * 
     * A.get_block(1, 1).to2dArray; // [[5, 6], [8, 9]]
     * A.get_block(0, 0, 1, 1).to2dArray; // [[1]]
     * A.get_block(1, 1, 2, 2).to2dArray; // [[5]]
     * A.get_block(0, 0, 2, 2).to2dArray; // [[1, 2], [4, 5]]
     */get_block(t,e,r=null,s=null){const[i,n]=this.shape;
/*if (!end_row)) {
            end_row = rows;
        }
            end_col = cols;
        }*/if(s=s??n,(r=r??i)<=t||s<=e)throw`\n                end_row must be greater than start_row, and \n                end_col must be greater than start_col, but\n                end_row = ${r}, start_row = ${t}, end_col = ${s}, and start_col = ${e}!`;const o=new Matrix(r-t,s-e,"zeros");for(let i=t,n=0;i<r;++i,++n)for(let t=e,r=0;t<s;++t,++r)o.set_entry(n,r,this.entry(i,t));return o;
//return new Matrix(end_row - start_row, end_col - start_col, (i, j) => this.entry(i + start_row, j + start_col));
}
/**
     * Returns a new array gathering entries defined by the indices given by argument.
     * @param {Array<Number>} row_indices - Array consists of indices of rows for gathering entries of this matrix 
     * @param {Array<Number>} col_indices  - Array consists of indices of cols for gathering entries of this matrix 
     * @returns {Matrix}
     */gather(t,e){const r=t.length,s=e.length,i=new Matrix(r,s);for(let s=0;s<r;++s){const n=t[s];for(let t=0;t<r;++t){const r=e[t];i.set_entry(s,t,this.entry(n,r))}}return i}
/**
     * Applies a function to each entry of the matrix.
     * @param {function} f function takes 2 parameters, the value of the actual entry and a value given by the function {@link v}. The result of {@link f} gets writen to the Matrix.
     * @param {function} v function takes 2 parameters for row and col, and returns a value witch should be applied to the colth entry of the rowth row of the matrix.
     */_apply_array(t,e){const r=this._data,[s,i]=this.shape;for(let n=0;n<s;++n){const s=n*i;for(let o=0;o<i;++o){const i=s+o;r[i]=t(r[i],e(n,o))}}return this}_apply_rowwise_array(t,e){return this._apply_array(e,((e,r)=>t[r]))}_apply_colwise_array(t,e){const r=this._data,[s,i]=this.shape;for(let n=0;n<s;++n){const s=n*i;for(let o=0;o<i;++o){const i=s+o;r[i]=e(r[i],t[n])}}return this}_apply(t,e){let r=this._data;if(t instanceof Matrix){let[s,i]=t.shape,[n,o]=this.shape;if(1===s){if(o!==i)throw"cols !== value_cols";for(let s=0;s<n;++s)for(let i=0;i<o;++i)r[s*o+i]=e(r[s*o+i],t.entry(0,i))}else if(1===i){if(n!==s)throw"rows !== value_rows";for(let s=0;s<n;++s)for(let i=0;i<o;++i)r[s*o+i]=e(r[s*o+i],t.entry(s,0))}else{if(n!=s||o!=i)throw"error";for(let s=0;s<n;++s)for(let i=0;i<o;++i)r[s*o+i]=e(r[s*o+i],t.entry(s,i))}}else if(Array.isArray(t)){let s=this._rows,i=this._cols;if(t.length===s)for(let n=0;n<s;++n)for(let s=0;s<i;++s)r[n*i+s]=e(r[n*i+s],t[n]);else{if(t.length!==i)throw"error";for(let n=0;n<s;++n)for(let s=0;s<i;++s)r[n*i+s]=e(r[n*i+s],t[s])}}else for(let s=0,i=this._rows*this._cols;s<i;++s)r[s]=e(r[s],t);return this}
/**
     * Clones the Matrix.
     * @returns {Matrix}
     */clone(){let t=new Matrix;return t._rows=this._rows,t._cols=this._cols,t._data=this._data.slice(0),t}mult(t){return this.clone()._apply(t,((t,e)=>t*e))}divide(t){return this.clone()._apply(t,((t,e)=>t/e))}add(t){return this.clone()._apply(t,((t,e)=>t+e))}sub(t){return this.clone()._apply(t,((t,e)=>t-e))}
/**
     * Returns the number of rows and columns of the Matrix.
     * @returns {Array} An Array in the form [rows, columns].
     */get shape(){return[this._rows,this._cols]}
/**
     * Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.
     * @param {Array} parameter - takes an Array in the form [rows, cols, value], where rows and cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row and col) which has to return a value for the colth entry of the rowth row.
     * @returns {Matrix}
     */set shape([t,e,r=(()=>0)]){this._rows=t,this._cols=e,this._data=new Float64Array(t*e);for(let s=0;s<t;++s)for(let t=0;t<e;++t)this._data[s*e+t]=r(s,t);return this}
/**
     * Returns the Matrix as a two-dimensional Array.
     * @returns {Array}
     */get to2dArray(){
/* const rows = this._rows;
        const cols = this._cols;
        let result = new Array(rows)
        for (let row = 0; row < rows; ++row) {
            let result_col = new Array(cols)
            for (let col = 0; col < cols; ++col) {
                result_col[col] = this.entry(row, col);
            }
            result[row] = result_col;
        }
        return result; */
return[...this.iterate_rows()]}
/**
     * Returns the diagonal of the Matrix.
     * @returns {Array}
     */get diag(){const t=this._rows,e=this._cols,r=Math.min(t,e);let s=new Float64Array(r);for(let t=0;t<r;++t)s[t]=this.entry(t,t);return s}
/**
     * Returns the mean of all entries of the Matrix.
     * @returns {float64}
     */get mean(){return this.sum/(this._rows*this._cols)}
/**
     * Returns the sum oof all entries of the Matrix.
     * @returns {number}
     */get sum(){return neumair_sum(this._data)}
/**
     * Returns the mean of each row of the matrix.
     * @returns {Array}
     */get meanRows(){const t=this._data,e=this._rows,r=this._cols;let s=Float64Array.from({length:e});for(let i=0;i<e;++i){s[i]=0;for(let e=0;e<r;++e)s[i]+=t[i*r+e];s[i]/=r}return s}
/** Returns the mean of each column of the matrix.
     * @returns {Array}
     */get meanCols(){const t=this._data,e=this._rows,r=this._cols;let s=Float64Array.from({length:r});for(let i=0;i<r;++i){s[i]=0;for(let n=0;n<e;++n)s[i]+=t[n*r+i];s[i]/=e}return s}static solve_CG(t,e,r,s=.001){const i=t.shape[0],n=e.shape[1];let o=new Matrix(i,0);for(let a=0;a<n;++a){const n=Matrix.from(e.col(a)).T;let h=new Matrix(i,1,(()=>r.random)),l=n.sub(t.dot(h)),_=l.clone();do{const e=t.dot(_),r=l.T.dot(l).entry(0,0)/_.T.dot(e).entry(0,0);h=h.add(_.mult(r));const s=l.sub(e.mult(r)),i=s.T.dot(s).entry(0,0)/l.T.dot(l).entry(0,0);_=s.add(_.mult(i)),l=s}while(Math.abs(l.mean)>s);o=o.concat(h,"horizontal")}return o}
/**
     * Solves the equation {@link A}x = {@link b}. Returns the result x.
     * @param {Matrix} A - Matrix or LU Decomposition
     * @param {Matrix} b - Matrix
     * @returns {Matrix}
     */static solve(t,e){let{L:r,U:s}="L"in t&&"U"in t?t:Matrix.LU(t),i=r.shape[0],n=e.clone();
// forward
for(let t=0;t<i;++t){for(let e=0;e<t-1;++e)n.set_entry(0,t,n.entry(0,t)-r.entry(t,e)*n.entry(1,e));n.set_entry(0,t,n.entry(0,t)/r.entry(t,t))}
// backward
for(let t=i-1;t>=0;--t){for(let e=i-1;e>t;--e)n.set_entry(0,t,n.entry(0,t)-s.entry(t,e)*n.entry(0,e));n.set_entry(0,t,n.entry(0,t)/s.entry(t,t))}return n}
/**
     * {@link L}{@link U} decomposition of the Matrix {@link A}. Creates two matrices, so that the dot product LU equals A.
     * @param {Matrix} A 
     * @returns {{L: Matrix, U: Matrix}} result - Returns the left triangle matrix {@link L} and the upper triangle matrix {@link U}.
     */static LU(t){const e=t.shape[0],r=new Matrix(e,e,"zeros"),s=new Matrix(e,e,"identity");for(let i=0;i<e;++i){for(let n=i;n<e;++n){let e=0;for(let t=0;t<i;++t)e+=r.entry(n,t)*s.entry(t,i);r.set_entry(n,i,t.entry(n,i)-e)}for(let n=i;n<e;++n){if(0===r.entry(i,i))return;let e=0;for(let t=0;t<i;++t)e+=r.entry(i,t)*s.entry(t,n);s.set_entry(i,n,(t.entry(i,n)-e)/r.entry(i,i))}}return{L:r,U:s}}
/**
     * Computes the {@link k} components of the SVD decomposition of the matrix {@link M}
     * @param {Matrix} M 
     * @param {int} [k=2] 
     * @returns {{U: Matrix, Sigma: Matrix, V: Matrix}}
     */static SVD(t,e=2){const r=t.T;let s=r.dot(t),i=t.dot(r),{eigenvectors:n,eigenvalues:o}=simultaneous_poweriteration(s,e),{eigenvectors:a}=simultaneous_poweriteration(i,e);return{U:a,Sigma:o.map((t=>Math.sqrt(t))),V:n};
//Algorithm 1a: Householder reduction to bidiagonal form:
/* const [m, n] = A.shape;
        let U = new Matrix(m, n, (i, j) => i == j ? 1 : 0);
        console.log(U.to2dArray)
        let V = new Matrix(n, m, (i, j) => i == j ? 1 : 0);
        console.log(V.to2dArray)
        let B = Matrix.bidiagonal(A.clone(), U, V);
        console.log(U,V,B)
        return { U: U, "Sigma": B, V: V }; */}}function distance_matrix(t,e=euclidean){let r=t.shape[0];
/* let D = new Array(n);
    for (let i = 0; i < n; ++i) {
        D[i] = new Float64Array(n);
    }
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            D[i][j] = D[j][i] = metric(A[i], A[j]);
        }
    } */const s=new Matrix(r,r);for(let i=0;i<r;++i){const n=t.row(i);for(let o=i+1;o<r;++o){const r=e(n,t.row(o));s.set_entry(i,o,r),s.set_entry(o,i,r)}}return s}function linspace(t,e,r=null){if(r||(r=Math.max(Math.round(e-t)+1,1)),r<2)return 1===r?[t]:[];let s=new Array(r);for(let i=r-=1;i>=0;--i)s[i]=(i*e+(r-i)*t)/r;return s}
//import { neumair_sum } from "../numerical/index";
function norm(t,e=euclidean){
//export default function(vector, p=2, metric = euclidean) {
let r=null;if(t instanceof Matrix){let[e,s]=t.shape;if(1===e)r=t.row(0);else{if(1!==s)throw"matrix must be 1d!";r=t.col(0)}}else r=t;let s=r.length,i=new Array(s);return i.fill(0),e(r,i);
/*let v;
    if (vector instanceof Matrix) {
        let [ rows, cols ] = v.shape;
        if (rows === 1) {
            v = vector.row(0);
        } else if (cols === 1) {
            v = vector.col(0);
        } else {
            throw "matrix must be 1d"
        }
    } else {
        v = vector;
    }
    return Math.pow(neumair_sum(v.map(e => Math.pow(e, p))), 1 / p)*/}class Randomizer{
// https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js
/**
     * Mersenne Twister
     * @param {*} _seed 
     */
constructor(t){return this._N=624,this._M=397,this._MATRIX_A=2567483615,this._UPPER_MASK=2147483648,this._LOWER_MASK=2147483647,this._mt=new Array(this._N),this._mti=this.N+1,this.seed=t||(new Date).getTime(),this}set seed(t){this._seed=t;let e=this._mt;for(e[0]=t>>>0,this._mti=1;this._mti<this._N;this._mti+=1){let t=this._mti,r=e[t-1]^e[t-1]>>>30;e[t]=(1812433253*((4294901760&r)>>>16)<<16)+1812433253*(65535&r)+t,e[t]>>>=0}}get seed(){return this._seed}get random(){return this.random_int*(1/4294967296)}get random_int(){let t,e=new Array(0,this._MATRIX_A);if(this._mti>=this._N){let r;this._mti==this._N+1&&(this.seed=5489);let s=this._N-this._M,i=this._M-this._N;for(r=0;r<s;++r)t=this._mt[r]&this._UPPER_MASK|this._mt[r+1]&this._LOWER_MASK,this._mt[r]=this._mt[r+this._M]^t>>>1^e[1&t];for(;r<this._N-1;++r)t=this._mt[r]&this._UPPER_MASK|this._mt[r+1]&this._LOWER_MASK,this._mt[r]=this._mt[r+i]^t>>>1^e[1&t];t=this._mt[this._N-1]&this._UPPER_MASK|this._mt[0]&this._LOWER_MASK,this._mt[this._N-1]=this._mt[this._M-1]^t>>>1^e[1&t],this._mti=0}return t=this._mt[this._mti+=1],t^=t>>>11,t^=t<<7&2636928640,t^=t<<15&4022730752,t^=t>>>18,t>>>0}choice(t,e){if(t instanceof Matrix){let[r,s]=t.shape;if(e>r)throw"n bigger than A!";let i=new Array(e),n=linspace(0,r-1);for(let t=0,r=n.length;t<e;++t,--r){let e=this.random_int%r;i[t]=n.splice(e,1)[0]}return i.map((e=>t.row(e)))}if(Array.isArray(t)||t instanceof Float64Array){let r=t.length;if(e>r)throw"n bigger than A!";let s=new Array(e),i=linspace(0,r-1);for(let t=0,r=i.length;t<e;++t,--r){let e=this.random_int%r;s[t]=i.splice(e,1)[0]}return s.map((e=>t[e]))}}static choice(t,e,r=19870307){let[s,i]=t.shape;if(e>s)throw"n bigger than A!";let n=new Randomizer(r),o=new Array(e),a=linspace(0,s-1);
/*let index_list = new Array(rows);
        for (let i = 0; i < rows; ++i) {
            index_list[i] = i;
        }*/
//let result = new Matrix(n, cols);
for(let t=0,r=a.length;t<e;++t,--r){let e=n.random_int%r;o[t]=a.splice(e,1)[0]}
//return result;
//return new Matrix(n, cols, (row, col) => A.entry(sample[row], col))
return o.map((e=>t.row(e)))}}function max(t){let e;for(const r of t)null!=r&&(e<r||void 0===e&&r>=r)&&(e=r);return e}
/**
 * @class
 * @alias Heap
 */class Heap{
/**
     * A heap is a datastructure holding its elements in a specific way, so that the top element would be the first entry of an ordered list.
     * @constructor
     * @memberof module:datastructure
     * @alias Heap
     * @param {Array=} elements - Contains the elements for the Heap. {@link elements} can be null.
     * @param {Function} [accessor = (d) => d] - Function returns the value of the element.
     * @param {("min"|"max"|Function)} [comparator = "min"] - Function returning true or false defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a Max_heap)
     * @returns {Heap}
     * @see {@link https://en.wikipedia.org/wiki/Binary_heap}
     */
constructor(t=null,e=(t=>t),r="min"){return t?Heap.heapify(t,e,r):(this._accessor=e,this._container=[],this._comparator="min"==r?(t,e)=>t<e:"max"==r?(t,e)=>t>e:r,this)}
/**
     * Creates a Heap from an Array
     * @param {Array|Set} elements - Contains the elements for the Heap.
     * @param {Function=} [accessor = (d) => d] - Function returns the value of the element.
     * @param {(String=|Function)} [comparator = "min"] - Function returning true or false defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a Max_heap)
     * @returns {Heap}
     */static heapify(t,e=(t=>t),r="min"){const s=new Heap(null,e,r),i=s._container;for(const r of t)i.push({element:r,value:e(r)});for(let e=Math.floor(t.length/2-1);e>=0;--e)s._heapify_down(e);return s}
/**
     * Swaps elements of container array.
     * @private
     * @param {Number} index_a 
     * @param {Number} index_b 
     */_swap(t,e){const r=this._container;[r[e],r[t]]=[r[t],r[e]]}
/**
     * @private
     */_heapify_up(){const t=this._container;let e=t.length-1;for(;e>0;){let r=Math.floor((e-1)/2);if(!this._comparator(t[e].value,t[r].value))break;this._swap(r,e),e=r}}
/**
     * Pushes the element to the heap.
     * @param {} element
     * @returns {Heap}
     */push(t){const e={element:t,value:this._accessor(t)};
//const node = new Node(element, value);
return this._container.push(e),this._heapify_up(),this}
/**
     * @private
     * @param {Number} [start_index = 0] 
     */_heapify_down(t=0){const e=this._container,r=this._comparator,s=e.length;let i=2*t+1,n=2*t+2,o=t;if(o>s)throw"index higher than length";i<s&&r(e[i].value,e[o].value)&&(o=i),n<s&&r(e[n].value,e[o].value)&&(o=n),o!==t&&(this._swap(t,o),this._heapify_down(o))}
/**
     * Removes and returns the top entry of the heap.
     * @returns {Object} Object consists of the element and its value (computed by {@link accessor}).
     */pop(){const t=this._container;if(0===t.length)return null;if(1===t.length)return t.pop();this._swap(0,t.length-1);const e=t.pop();return this._heapify_down(),e}
/**
     * Returns the top entry of the heap without removing it.
     * @returns {Object} Object consists of the element and its value (computed by {@link accessor}).
     */get first(){return this._container.length>0?this._container[0]:null}
/**
     * Yields the raw data
     * @yields {Object} Object consists of the element and its value (computed by {@link accessor}).
     */*iterate(){for(let t=0,e=this._container.length;t<e;++t)yield this._container[t].element}
/**
     * Returns the heap as ordered array.
     * @returns {Array} Array consisting the elements ordered by {@link comparator}.
     */toArray(){return this.data().sort(((t,e)=>this._comparator(t,e)?-1:0))}
/**
     * Returns elements of container array.
     * @returns {Array} Array consisting the elements.
     */data(){return this._container.map((t=>t.element))}
/**
     * Returns the container array.
     * @returns {Array} The container array.
     */raw_data(){return this._container}
/**
     * The size of the heap.
     * @returns {Number}
     */get length(){return this._container.length}
/**
     * Returns false if the the heap has entries, true if the heap has no entries.
     * @returns {Boolean}
     */get empty(){return 0===this.length}}
/**
 * @class
 * @alias DisjointSet
 * @see {@link https://en.wikipedia.org/wiki/Disjoint-set_data_structure}
 */class DisjointSet{
/**
     * @constructor
     * @alias DisjointSet
     * @memberof module:datastructure
     * @param {Array=} elements 
     * @returns {DisjointSet}
     */
constructor(t=null){if(this._list=new Set,t)for(const e of t)this.make_set(e);return this}make_set(t){const e=this._list;return e.has(t)||(e.add(t),t.__disjoint_set={},t.__disjoint_set.parent=t,t.__disjoint_set.children=new Set([t]),t.__disjoint_set.size=1),this}find(t){return this._list.has(t)?t.__disjoint_set.parent!==t?(t.__disjoint_set.children.add(...t),t.__disjoint_set.parent=this.find(t.__disjoint_set.parent),t.__disjoint_set.parent):t:null}union(t,e){let r=this.find(t),s=this.find(e);return r===s||(r.__disjoint_set.size<s.__disjoint_set.size&&([r,s]=[s,r]),s.__disjoint_set.parent=r,
// keep track of children?
s.__disjoint_set.children.forEach(r.__disjoint_set.children.add,r.__disjoint_set.children),r.__disjoint_set.size+=s.__disjoint_set.size),this}}
/**
 * @class
 * @alias BallTree
 */class BallTree{
/**
     * Generates a BallTree with given {@link elements}.
     * @constructor
     * @memberof module:knn
     * @alias BallTree
     * @param {Array=} elements - Elements which should be added to the BallTree
     * @param {Function} [metric = euclidean] metric to use: (a, b) => distance
     * @see {@link https://en.wikipedia.org/wiki/Ball_tree}
     * @see {@link https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js}
     * @returns {BallTree}
     */
constructor(t=null,e=euclidean){return this._Node=class{constructor(t,e=null,r=null,s=null){this.pivot=t,this.child1=e,this.child2=r,this.radius=s}},this._Leaf=class{constructor(t){this.points=t}},this._metric=e,t&&this.add(t),this}
/**
     * 
     * @param {Array<*>} elements - new elements.
     * @returns {BallTree}
     */add(t){return t=t.map(((t,e)=>({index:e,element:t}))),this._root=this._construct(t),this}
/**
     * @private
     * @param {Array<*>} elements 
     * @returns {Node} root of balltree.
     */_construct(t){if(1===t.length)return new this._Leaf(t);{let e,r=this._greatest_spread(t),s=t.sort(((t,e)=>t.element[r]-e.element[r])),i=s.length,n=Math.floor(i/2),o=t[n],a=s.slice(0,n),h=s.slice(n,i),l=Math.max(...t.map((t=>this._metric(o.element,t.element))));return e=a.length>0&&h.length>0?new this._Node(o,this._construct(a),this._construct(h),l):new this._Leaf(t),e}}
/**
     * @private
     * @param {Node} B 
     * @returns {Number}
     */_greatest_spread(t){let e=t[0].element.length,r=new Array(e);for(let t=0;t<e;++t)r[t]=[1/0,-1/0];let s=t.reduce(((t,r)=>{for(let s=0;s<e;++s)t[s][0]=Math.min(t[s][0],r.element[s]),t[s][1]=Math.max(t[s][1],r.element[s]);return t}),r);s=s.map((t=>t[1]-t[0]));let i=0;for(let t=0;t<e;++t)i=s[t]>s[i]?t:i;return i}
/**
     * 
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @returns {Heap} - Heap consists of the {@link k} nearest neighbors.
     */search(t,e=5){return this._search(t,e,new Heap(null,(e=>this._metric(e.element,t)),"max"),this._root)}
/**
     * @private
     * @param {*} t - query element.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @param {Heap} Q - Heap consists of the currently found {@link k} nearest neighbors.
     * @param {Node|Leaf} B 
     */_search(t,e,r,s){
// B is Node
if(r.length>=e&&s.pivot&&s.radius&&this._metric(t,s.pivot.element)-s.radius>=r.first.value)return r;
// B is leaf
if(s.child1&&this._search(t,e,r,s.child1),s.child2&&this._search(t,e,r,s.child2),s.points)for(let t=0,i=s.points.length;t<i;++t){let i=s.points[t];e>r.length?r.push(i):(r.push(i),r.pop())}return r}}
/**
 * @class
 * @alias KNN
 */class KNN{
/**
     * Generates a KNN list with given {@link elements}.
     * @constructor
     * @memberof module:knn
     * @alias KNN
     * @param {Array=} elements - Elements which should be added to the KNN list
     * @param {Function|"precomputed"} [metric = euclidean] metric is either precomputed or a function to use: (a, b) => distance
     * @returns {KNN}
     */
constructor(t=null,e=euclidean){this._metric=e,this._elements=t instanceof Matrix?t:Matrix.from(t);const r=this._elements.shape[0];this._D="precomputed"===e?this._elements.clone():distance_matrix(this._elements,e),this.KNN=[];for(let t=0;t<r;++t){const e=this._D.row(t),s=new Heap(null,(t=>t.value),"min");for(let t=0;t<r;++t)s.push({value:e[t],index:t});this.KNN.push(s)}}
/**
     * 
     * @param {Array|Number} t - query element or index.
     * @param {Number} [k = 5] - number of nearest neighbors to return.
     * @returns {Heap} - Heap consists of the {@link k} nearest neighbors.
     */search(t,e=5){const r=this._metric,s=this.KNN;let i;if(Array.isArray(t)){if("precomputed"==this._metric)throw"Search by query element is only possible when not using a precomputed distance matrix!";const e=this._elements,n=s.length;let o=null,a=1/0;for(let s=0;s<n;++s){const i=r(t,e.row(s));i<a&&(o=s,a=i)}i=s[o]}else Number.isInteger(t)&&(i=s[t]);let n=[];for(let t=0;t<e;++t)n.push(i.pop());return n.forEach((t=>i.push(t.element))),n}}
/**
 * Computes the QR Decomposition of the Matrix {@link A} using Gram-Schmidt process.
 * @memberof module:linear_algebra
 * @alias qr
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process}
 */function qr(t){const[e,r]=t.shape,s=new Matrix(e,r,"identity"),i=new Matrix(r,r,0);for(let n=0;n<r;++n){let r=t.col(n);for(let t=0;t<n;++t){const e=s.col(t),o=neumair_sum(e.map(((t,e)=>t*r[e])));i.set_entry(t,n,o),r=r.map(((t,r)=>t-o*e[r]))}const o=norm(r,euclidean);for(let t=0;t<e;++t)s.set_entry(t,n,r[t]/o);i.set_entry(n,n,o)}return{R:i,Q:s}}
/**
 * Computes the {@link k} biggest Eigenvectors and Eigenvalues from Matrix {@link A} with the QR-Algorithm.
 * @param {Matrix} A - The Matrix
 * @param {Number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {Number} [max_iterations=100] - The number of maxiumum iterations the algorithm should run.
 * @param {Number|Randomizer} [seed=1212] - The seed value or a randomizer used in the algorithm.
 * @returns {{eigenvalues: Array, eigenvectors: Array}} - The {@link k} biggest eigenvectors and eigenvalues of Matrix {@link A}.
 */function simultaneous_poweriteration$1(t,e=2,r=100,s=1212){const i=s instanceof Randomizer?s:new Randomizer(s);t instanceof Matrix||(t=Matrix.from(t));const n=t.shape[0];let{Q:o,R:a}=qr(new Matrix(n,e,(()=>i.random)));for(;r--;){const e=a.clone(),s=qr(t.dot(o));o=s.Q,a=s.R,neumair_sum(a.sub(e).diag)/n<1e-12&&(r=0)}return{eigenvalues:a.diag,eigenvectors:o.transpose().to2dArray}}
/**
 * @class
 * @alias DR
 */class DR$1{
//static parameter_list = [];
get parameter_list(){return this._parameter_list}set parameter_list(t){return this._parameter_list=t,this}
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias DR
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data. 
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {seed} [seed=1987] - the seed value for the random number generator.
     * @returns {DR}
     */constructor(t,e=2,r=euclidean,s=1212){if(Array.isArray(t))this._type="array",this.X=Matrix.from(t);else{if(!(t instanceof Matrix))throw"no valid type for X";this._type="matrix",this.X=t}return[this._N,this._D]=this.X.shape,this._d=e,this._metric=r,this._seed=s,this._randomizer=new Randomizer(s),this._is_initialized=!1,this}
/**
     * Set and get parameters
     * @param {String} name - name of the parameter.
     * @param {Number} [value = null] - value of the parameter to set, if null then return actual parameter value.
     */parameter(t,e=null){if(-1===this.parameter_list.findIndex((e=>e===t)))throw t+" is not a valid parameter!";return e?(this["_"+t]=e,this):this["_"+t]}
/**
     * Alias for 'parameter'.
     * @param {String} name 
     * @param {Number} value 
     */para(t,e=null){return this.parameter(t,e)}
/**
     * Alias for 'parameter'.
     * @param {String} name 
     * @param {Number} value 
     */p(t,e=null){return this.parameter(t,e)}
/**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */transform(){return this.check_init(),this.Y}*generator(){return this.transform()}check_init(){this._is_initialized||"function"!=typeof this.init||(this.init(),this._is_initialized=!0)}
/**
     * @returns {Matrix} Returns the projection.
     */get projection(){return"matrix"===this._type?this.Y:this.Y.to2dArray}async transform_async(){return this.transform()}static transform(...t){return new this(...t).transform()}static async transform_async(...t){return this.transform(...t)}static*generator(...t){const e=new this(...t).generator();for(const t of e)yield t}}
/**
 * @class
 * @alias PCA
 * @augments DR
 */class PCA extends DR$1{
/**
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias PCA 
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @returns {PCA}
     */
constructor(t,e=2){return super(t,e),this}
/**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */transform(){let t=this.X,e=t.shape[1],r=new Matrix(e,e,"center"),s=t.dot(r),i=s.transpose().dot(s),{eigenvectors:n}=simultaneous_poweriteration$1(i,this._d);return n=Matrix.from(n).transpose(),this.Y=t.dot(n),this.projection}}
/**
 * @class
 * @alias MDS
 */class MDS extends DR$1{
/**
     * Classical MDS.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias MDS
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
constructor(t,e=2,r=euclidean,s=1212){return super(t,e,r,s),this}
/**
     * Transforms the inputdata {@link X} to dimensionality {@link d}.
     * @returns {Matrix|Array}
     */transform(){const t=this.X,e=t.shape[0],r=this._metric,s="precomputed"===r?t:distance_matrix(t,r),i=s.meanCols,n=s.meanRows,o=s.mean;this._d_X=s;const a=new Matrix(e,e,((t,e)=>s.entry(t,e)-i[t]-n[e]+o)),{eigenvectors:h}=simultaneous_poweriteration$1(a,this._d);return this.Y=Matrix.from(h).transpose(),this.projection}
/**
     * @returns {Number} - the stress of the projection.
     */stress(){const t=this.X.shape[0],e=this.Y,r=this._d_X,s=new Matrix;s.shape=[t,t,(t,r)=>t<r?euclidean(e.row(t),e.row(r)):s.entry(r,t)];let i=0,n=0;for(let e=0;e<t;++e)for(let o=e+1;o<t;++o)i+=Math.pow(r.entry(e,o)-s.entry(e,o),2),n+=Math.pow(r.entry(e,o),2);return Math.sqrt(i/n)}}
/**
 * @class
 * @alias ISOMAP
 */class ISOMAP extends DR$1{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias ISOMAP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} neighbors - the number of neighbors {@link ISOMAP} should use to project the data.
     * @param {Number} [d = 2] - the dimensionality of the projection. 
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points. 
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
constructor(t,e,r=2,s=euclidean,i=1212){return super(t,r,s,i),super.parameter_list=["k"],this.parameter("k",Math.min(e??Math.max(Math.floor(this.X.shape[0]/10),2),this._N-1)),this}
/**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */transform(){this.check_init();const t=this.X,e=this._N,r=this._metric,s=new Matrix;s.shape=[e,e,(e,i)=>e<=i?r(t.row(e),t.row(i)):s.entry(i,e)];const i=[];for(let t=0;t<e;++t){const r=[];for(let i=0;i<e;++i)r.push({index:i,distance:s.entry(t,i)});const n=new Heap(r,(t=>t.distance),"min");i.push(n.toArray().slice(1,this._k+1))}
/*D = dijkstra(kNearestNeighbors);*/
// compute shortest paths
// TODO: make extern
/** @see {@link https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm} */const n=new Matrix(e,e,((t,e)=>{const r=i[t].find((t=>t.index===e));return r?r.distance:1/0}));for(let t=0;t<e;++t)for(let r=0;r<e;++r)for(let s=0;s<e;++s)n.set_entry(t,r,Math.min(n.entry(t,r),n.entry(t,s)+n.entry(s,r)));let o=new Float64Array(e),a=new Float64Array(e),h=0,l=new Matrix(e,e,((t,e)=>{let r=n.entry(t,e);return r=r===1/0?0:r,o[t]+=r,a[e]+=r,h+=r,r}));o=o.map((t=>t/e)),a=a.map((t=>t/e)),h/=e**2;const _=new Matrix(e,e,((t,e)=>l.entry(t,e)-o[t]-a[e]+h)),{eigenvectors:c}=simultaneous_poweriteration$1(_,this._d);
// compute d eigenvectors
// return embedding
return this.Y=Matrix.from(c).transpose(),this.projection}}
/**
 * @class
 * @alias FASTMAP
 */class FASTMAP extends DR$1{
/**
     * FastMap: a fast algorithm for indexing, data-mining and visualization of traditional and multimedia datasets
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias FASTMAP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {FASTMAP}
     * @see {@link https://doi.org/10.1145/223784.223812}
     */
constructor(t,e=2,r=euclidean,s=1212){return super(t,e,r,s),this}
/**
     * Chooses two points which are the most distant in the actual projection.
     * @private
     * @param {function} dist 
     * @returns {Array} An array consisting of first index, second index, and distance between the two points.
     */_choose_distant_objects(t){const e=this.X.shape[0];let r=this._randomizer.random_int%e-1,s=null,i=-1/0;for(let n=0;n<e;++n){const e=t(r,n);e>i&&(i=e,s=n)}i=-1/0;for(let n=0;n<e;++n){const e=t(s,n);e>i&&(i=e,r=n)}return[r,s,i]}
/**
     * Computes the projection.
     * @returns {Matrix} The {@link d}-dimensional projection of the data matrix {@link X}.
     */transform(){const t=this.X,e=t.shape[0],r=this._d,s=this._metric,i=new Matrix(e,r,0);let dist=(e,r)=>s(t.row(e),t.row(r));for(let t=0;t<r;++t){let r=dist;
// choose pivot objects
const[s,n,o]=this._choose_distant_objects(dist);
// record id of pivot objects
//PA[0].push(a_index);
//PA[1].push(b_index);
/* if (d_ab === 0) {
                // because all inter-object distances are zeros
                for (let i = 0; i < N; ++i) {
                    Y.set_entry(i, _col, 0);
                }
            } else { */if(0!==o){
// project the objects on the line (O_a, O_b)
for(let r=0;r<e;++r){const e=(dist(s,r)**2+o**2-dist(n,r)**2)/(2*o);i.set_entry(r,t,e)}
// consider the projections of the objects on a
// hyperplane perpendicluar to the line (a, b);
// the distance function D'() between two 
// projections is given by Eq.4
dist=(e,s)=>Math.sqrt(r(e,s)**2-(i.entry(e,t)-i.entry(s,t))**2)}}
// return embedding
return this.Y=i,this.projection}}
/**
 * @class
 * @alias LDA
 */class LDA extends DR$1{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LDA
     * @param {Matrix} X - the high-dimensional data.
     * @param {Array} labels - the label / class of each data point.
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
constructor(t,e,r=2,s=euclidean,i=1212){return super(t,r,s,i),super.parameter_list=["labels"],this.parameter("labels",e),this}
/**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */transform(){let t=this.X,[e,r]=t.shape,s=this._labels,i={},n=0;s.forEach(((e,r)=>{e in i?(i[e].count++,i[e].rows.push(t.row(r))):i[e]={id:n++,count:1,rows:[t.row(r)]}}));
// create X_mean and vector means;
let o=t.mean,a=new Matrix(n,r);for(let t in i){let e=Matrix.from(i[t].rows).meanCols;for(let s=0;s<r;++s)a.set_entry(i[t].id,s,e[s])}
// scatter_between
let h=new Matrix(r,r);for(let t in i){let e=a.row(i[t].id),s=new Matrix(r,1,(t=>e[t]-o)),n=i[t].count;h=h.add(s.dot(s.transpose()).mult(n))}
// scatter_within
let l=new Matrix(r,r);for(let t in i){let e=a.row(i[t].id),s=new Matrix(r,1,(t=>e[t])),n=i[t].rows;for(let e=0,o=i[t].count;e<o;++e){let t=new Matrix(r,1,((t,r)=>n[e][t]-s.entry(t,0)));l=l.add(t.dot(t.transpose()))}}let{eigenvectors:_}=simultaneous_poweriteration$1(l.inverse().dot(h),this._d);
// return embedding
return _=Matrix.from(_).transpose(),this.Y=t.dot(_),this.projection}}
/**
 * @class
 * @alias LLE
 */class LLE extends DR$1{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LLE
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     */
constructor(t,e,r=2,s=euclidean,i=1212){return super(t,r,s,i),super.parameter_list=["k"],this.parameter("k",Math.min(e??Math.max(Math.floor(this._N/10),2),this._N-1)),this}
/**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */transform(){const t=this.X,e=this._d,r=this._N,s=this._D,i=this.parameter("k"),n=k_nearest_neighbors(t,i,null,this._metric),o=new Matrix(i,1,1),a=new Matrix(r,r);for(let e=0;e<r;++e){const r=n[e],h=new Matrix(i,s,((s,i)=>t.entry(r[s].j,i)-t.entry(e,i))),l=h.dot(h.T);if(i>s){const t=neumair_sum(l.diag)/1e3;for(let e=0;e<i;++e)l.set_entry(e,e,l.entry(e,e)+t)}
// reconstruct;
let _=Matrix.solve_CG(l,o,this._randomizer);_=_.divide(_.sum);for(let t=0;t<i;++t)a.set_entry(e,r[t].j,_.entry(t,0))}
// comp embedding
const h=new Matrix(r,r,"identity").sub(a),l=h.T.dot(h),{eigenvectors:_}=simultaneous_poweriteration$1(l.T.inverse(),e+1);
// return embedding
return this.Y=Matrix.from(_.slice(1,1+e)).T,this.projection}}
/**
 * @class
 * @alias LTSA
 */class LTSA extends DR$1{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LTSA
     * @param {Matrix} X - the high-dimensional data.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
     */
constructor(t,e,r=2,s=euclidean,i=1212){if(super(t,r,s,i),super.parameter_list=["k"],this.parameter("k",Math.min(e??Math.max(Math.floor(this._N/10),2),this._N-1)),this._D<=r)throw`Dimensionality of X (D = ${this._D}) must be greater than the required dimensionality of the result (d = ${r})!`;return this}
/**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */transform(){const t=this.X,e=this._d,[r,s]=t.shape,i=this.parameter("k"),n=k_nearest_neighbors(t,i,null,this._metric),o=new Matrix(s,s,"center"),a=new Matrix(r,r,0);for(let s=0;s<r;++s){
// 1.2 compute the d largest eigenvectors of the correlation matrix
const r=[s,...n[s].map((t=>t.j))];let h=Matrix.from(r.map((e=>t.row(e))));
// center X_i
h=h.dot(o);
// correlation matrix
const l=h.dot(h.transpose()),{eigenvectors:_}=simultaneous_poweriteration$1(l,e),c=Matrix.from(_),u=c.transpose().dot(c).add(1/Math.sqrt(i+1));for(let t=0;t<i+1;++t)for(let e=0;e<i+1;++e)a.set_entry(r[t],r[e],a.entry(r[t],r[e])-(t===e?1:0)+u.entry(t,e))}
// 3. Aligning global coordinates
const{eigenvectors:h}=simultaneous_poweriteration$1(a,e+1);
// return embedding
return this.Y=Matrix.from(h.slice(1)).transpose(),this.projection}}
/**
 * @class
 * @alias TSNE
 */class TSNE extends DR$1{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TSNE
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [perplexity = 50] - perplexity.
     * @param {Number} [epsilon = 10] - learning parameter.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {TSNE}
     */
constructor(t,e=50,r=10,s=2,i=euclidean,n=1212){return super(t,s,i,n),super.parameter_list=["perplexity","epsilon"],[this._N,this._D]=this.X.shape,this.parameter("perplexity",Math.min(e,this._N-1)),this.parameter("epsilon",r),this._iter=0,this.Y=new Matrix(this._N,this._d,(()=>this._randomizer.random)),this}
/**
     * 
     * @param {Matrix} distance_matrix - accepts a precomputed distance matrix
     * @returns {TSNE}
     */init(t=null){
// init
const e=Math.log(this._perplexity),r=this._N,s=this._D,i=this._metric,n=this.X;let o;if(t)o=t;else{o=new Matrix(r,r);for(let t=0;t<r;++t){const e=n.row(t);for(let s=t+1;s<r;++s){const r=i(e,n.row(s));o.set_entry(t,s,r),o.set_entry(s,t,r)}}}const a=new Matrix(r,r,"zeros");this._ystep=new Matrix(r,s,"zeros"),this._gains=new Matrix(r,s,1);
// search for fitting sigma
let h=new Array(r).fill(0);for(let t=0;t<r;++t){let s=-1/0,i=1/0,n=1,l=!1,_=0;for(;!l;){let a=0;for(let e=0;e<r;++e){let r=Math.exp(-o.entry(t,e)*n);t===e&&(r=0),h[e]=r,a+=r}let c=0;for(let t=0;t<r;++t){let e=0===a?0:h[t]/a;h[t]=e,e>1e-7&&(c-=e*Math.log(e))}c>e?(s=n,n=i===1/0?2*n:(n+i)/2):(i=n,n=s===-1/0?n/2:(n+s)/2),++_,Math.abs(c-e)<1e-4&&(l=!0),_>=50&&(l=!0)}for(let e=0;e<r;++e)a.set_entry(t,e,h[e])}
//compute probabilities
const l=new Matrix(r,r,"zeros"),_=2*r;for(let t=0;t<r;++t)for(let e=t;e<r;++e){const r=Math.max((a.entry(t,e)+a.entry(e,t))/_,1e-100);l.set_entry(t,e,r),l.set_entry(e,t,r)}return this._P=l,this}
/**
     * 
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Array<Array>} - the projection.
     */transform(t=500){this.check_init();for(let e=0;e<t;++e)this.next();return this.projection}
/**
     * 
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Array<Array>} - the projection.
     */*generator(t=500){this.check_init();for(let e=0;e<t;++e)this.next(),yield this.projection;return this.projection}
/**
     * performs a optimization step
     * @private
     * @returns {Matrix}
     */next(){const t=++this._iter,e=this._P,r=this._ystep,s=this._gains,i=this._N,n=this._epsilon,o=this._d;let a=this.Y;
//calc cost gradient;
const h=t<100?4:1,l=new Matrix(i,i,"zeros");
// compute Q dist (unnormalized)
let _=0;for(let t=0;t<i;++t)for(let e=t+1;e<i;++e){let r=0;for(let s=0;s<o;++s){const i=a.entry(t,s)-a.entry(e,s);r+=i*i}const s=1/(1+r);l.set_entry(t,e,s),l.set_entry(e,t,s),_+=2*s}
// normalize Q dist
const c=new Matrix(i,i,0);for(let t=0;t<i;++t)for(let e=t+1;e<i;++e){const r=Math.max(l.entry(t,e)/_,1e-100);c.set_entry(t,e,r),c.set_entry(e,t,r)}const u=new Matrix(i,o,"zeros");for(let t=0;t<i;++t)for(let r=0;r<i;++r){const s=4*(h*e.entry(t,r)-c.entry(t,r))*l.entry(t,r);for(let e=0;e<o;++e)u.set_entry(t,e,u.entry(t,e)+s*(a.entry(t,e)-a.entry(r,e)))}
// perform gradient step
let d=new Float64Array(o);for(let e=0;e<i;++e)for(let i=0;i<o;++i){const o=u.entry(e,i),h=r.entry(e,i),l=s.entry(e,i);let _=Math.sign(o)===Math.sign(h)?.8*l:l+.2;_<.01&&(_=.01),s.set_entry(e,i,_);const c=(t<250?.5:.8)*h-n*_*o;r.set_entry(e,i,c),a.set_entry(e,i,a.entry(e,i)+c),d[i]+=a.entry(e,i)}for(let t=0;t<i;++t)for(let e=0;e<2;++e)a.set_entry(t,e,a.entry(t,e)-d[e]/i);return this.Y}}
// http://optimization-js.github.io/optimization-js/optimization.js.html#line438
function powell(t,e,r=300){const s=e.length;let i=.001,n=1e4,o=e.slice(),a=t(o),h=!1;for(;r-- >=0&&!h;){h=!0;for(let e=0;e<s;++e){o[e]+=1e-6;let r=t(o);o[e]-=1e-6;let s=(r-a)/1e-6;Math.abs(s)>.01&&(h=!1),o[e]-=i*s,a=t(o)}i*=n>=a?1.05:.4,n=a}return o}
/**
 * @class
 * @alias UMAP
 */class UMAP extends DR$1{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias UMAP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [n_neighbors = 15] - size of the local neighborhood.
     * @param {Number} [local_connectivity = 1] - number of nearest neighbors connected in the local neighborhood.
     * @param {Number} [min_dist = 1] - controls how tightly points get packed together.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points in the high-dimensional space.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {UMAP}
     */
constructor(t,e=15,r=1,s=1,i=2,n=euclidean,o=1212){return super(t,i,n,o),super.parameter_list=["n_neighbors","local_connectivity","min_dist"],[this._N,this._D]=this.X.shape,e=Math.min(this._N-1,e),this.parameter("n_neighbors",e),this.parameter("local_connectivity",Math.min(r,e-1)),this.parameter("min_dist",s),this._iter=0,this._spread=1,this._set_op_mix_ratio=1,this._repulsion_strength=1,this._negative_sample_rate=5,this._n_epochs=350,this._initial_alpha=1,this.Y=new Matrix(this._N,this._d,(()=>this._randomizer.random)),this}
/**
     * @private
     * @param {Number} spread 
     * @param {Number} min_dist 
     * @returns {Array}
     */_find_ab_params(t,e){const r=linspace(0,3*t,300),s=linspace(0,3*t,300);for(let i=0,n=r.length;i<n;++i){const n=r[i];s[i]=n<e?1:Math.exp(-(n-e)/t)}return powell((t=>{const e=linspace(1,300).map(((e,i)=>{return s[i]-(n=r[i],o=t[0],a=t[1],1/(1+o*Math.pow(n,2*a)));var n,o,a}));return Math.sqrt(neumair_sum(e.map((t=>t*t))))}),[1,1])}
/**
     * @private
     * @param {Array<Array>} distances 
     * @param {Array<Number>} sigmas 
     * @param {Array<Number>} rhos 
     * @returns {Array}
     */_compute_membership_strengths(t,e,r){for(let s=0,i=t.length;s<i;++s)for(let i=0,n=t[s].length;i<n;++i){const n=t[s][i].value-r[s];t[s][i].value=n>0?Math.exp(-n/e[s]):1}return t}
/**
     * @private
     * @param {KNN|BallTree} knn 
     * @param {Number} k 
     * @returns {Object}
     */_smooth_knn_dist(t,e){const r=1e-5,s=.001,i=this._local_connectivity,n=Math.log2(e),o=[],a=[],h=this.X,l=h.shape[0],_=[];if("precomputed"===this._metric)for(let r=0;r<l;++r)_.push(t.search(r,e).reverse());else for(const r of h)_.push(t.search(r,e).raw_data().reverse());for(let t=0;t<l;++t){let h=0,l=1/0,c=1;const u=_[t],d=u.filter((t=>t.value>0)),f=d.length;if(f>=i){const e=Math.floor(i),s=i-e;e>0?(o.push(d[e-1]),s>r&&(o[t].value+=s*(d[e].value-d[e-1]))):o[t].value=s*d[0].value}else f>0&&(o[t]=d[f-1].value);for(let s=0;s<64;++s){let s=0;for(let r=0;r<e;++r){const e=u[r].value-o[t];s+=e>0?Math.exp(-e/c):1}if(Math.abs(s-n)<r)break;s>n?[l,c]=[c,(h+l)/2]:[h,c]=l===1/0?[c,2*c]:[c,(h+l)/2]}a[t]=c;const p=u.reduce(((t,e)=>t+e.value),0)/u.length;
//let mean_d = null;
if(o[t]>0)a[t]<s*p&&(a[t]=s*p);else{const e=_.reduce(((t,e)=>t+e.reduce(((t,e)=>t+e.value),0)/e.length));a[t]>s*e&&(a[t]=s*e)}}return{distances:_,sigmas:a,rhos:o}}
/**
     * @private
     * @param {Matrix} X 
     * @param {Number} n_neighbors 
     * @returns {Matrix}
     */_fuzzy_simplicial_set(t,e){const r=t.shape[0],s=this._metric,i="precomputed"===s?new KNN(t,"precomputed"):new BallTree(t.to2dArray,s);let{distances:n,sigmas:o,rhos:a}=this._smooth_knn_dist(i,e);n=this._compute_membership_strengths(n,o,a);const h=new Matrix(r,r,"zeros");for(let t=0;t<r;++t){const e=n[t];for(let r=0;r<e.length;++r)h.set_entry(t,e[r].element.index,e[r].value)}const l=h.T,_=h.mult(l);return h.add(l).sub(_).mult(this._set_op_mix_ratio).add(_.mult(1-this._set_op_mix_ratio))}
/**
     * @private
     * @param {Number} n_epochs 
     * @returns {Array}
     */_make_epochs_per_sample(t){const e=this._weights,r=new Float32Array(e.length).fill(-1),s=max(e),i=e.map((e=>t*(e/s)));for(let e=0;e<r.length;++e)i[e]>0&&(r[e]=Math.round(t/i[e]));return r}
/**
     * @private
     * @param {Matrix} graph 
     * @returns {Object}
     */_tocoo(t){const e=[],r=[],s=[],[i,n]=t.shape;for(let o=0;o<i;++o)for(let i=0;i<n;++i){const n=t.entry(o,i);0!==n&&(e.push(o),r.push(i),s.push(n))}return{rows:e,cols:r,data:s}}
/**
     * Computes all necessary 
     * @returns {UMAP}
     */init(){const[t,e]=this._find_ab_params(this._spread,this._min_dist);this._a=t,this._b=e,this._graph=this._fuzzy_simplicial_set(this.X,this._n_neighbors);const{rows:r,cols:s,data:i}=this._tocoo(this._graph);return this._head=r,this._tail=s,this._weights=i,this._epochs_per_sample=this._make_epochs_per_sample(this._n_epochs),this._epochs_per_negative_sample=this._epochs_per_sample.map((t=>t*this._negative_sample_rate)),this._epoch_of_next_sample=this._epochs_per_sample.slice(),this._epoch_of_next_negative_sample=this._epochs_per_negative_sample.slice(),this}set local_connectivity(t){this._local_connectivity=t}get local_connectivity(){return this._local_connectivity}set min_dist(t){this._min_dist=t}get min_dist(){return this._min_dist}graph(){return this.check_init(),{cols:this._head,rows:this._tail,weights:this._weights}}
/**
     * 
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */transform(t=350){this._n_epochs!=t&&(this._n_epochs=t,this.init()),this.check_init();for(let e=0;e<t;++e)this.next();return this.projection}
/**
     * 
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */*generator(t=350){this._n_epochs!=t&&(this._n_epochs=t,this.init()),this.check_init();for(let e=0;e<t;++e)this.next(),yield this.projection;return this.projection}
/**
     * @private
     * @param {Number} x 
     * @returns {Number}
     */_clip(t){return t>4?4:t<-4?-4:t}
/**
     * performs the optimization step.
     * @private
     * @param {Matrix} head_embedding 
     * @param {Matrix} tail_embedding 
     * @param {Matrix} head 
     * @param {Matrix} tail 
     * @returns {Matrix}
     */_optimize_layout(t,e,r,s){const{_d:i,_alpha:n,_repulsion_strength:o,_a:a,_b:h,_epochs_per_sample:l,_epochs_per_negative_sample:_,_epoch_of_next_negative_sample:c,_epoch_of_next_sample:u,_clip:d}=this,f=s.length;for(let p=0,m=l.length;p<m;++p)if(u[p]<=this._iter){const m=r[p],y=s[p],w=t.row(m),g=e.row(y),x=euclidean_squared(w,g);let M=0;x>0&&(M=-2*a*h*Math.pow(x,h-1)/(a*Math.pow(x,h)+1));for(let r=0;r<i;++r){const s=d(M*(w[r]-g[r]))*n,i=w[r]+s,o=g[r]-s;w[r]=i,g[r]=o,t.set_entry(m,r,i),e.set_entry(y,r,o)}u[p]+=l[p];const A=(this._iter-c[p])/_[p];for(let r=0;r<A;++r){const r=Math.floor(this._randomizer.random*f),l=e.row(s[r]),_=euclidean_squared(w,l);let c=0;if(_>0)c=2*o*h/((.01+_)*(a*Math.pow(_,h)+1));else if(m===r)continue;for(let o=0;o<i;++o){const i=d(c*(w[o]-l[o]))*n,a=w[o]+i,h=l[o]-i;w[o]=a,l[o]=h,t.set_entry(m,o,a),e.set_entry(s[r],o,h)}}c[p]+=A*_[p]}return t}
/**
     * @private
     * @returns {Matrix}
     */next(){let t=++this._iter,e=this.Y;return this._alpha=this._initial_alpha*(1-t/this._n_epochs),this.Y=this._optimize_layout(e,e,this._head,this._tail),this.Y}}
/**
 * @class
 * @alias TriMap
 */class TriMap extends DR$1{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TriMap
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [weight_adj = 500] - scaling factor.
     * @param {Number} [c = 5] - number of triplets multiplier.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {TriMap}
     * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
     * @see {@link https://github.com/eamid/trimap}
     */
constructor(t,e=500,r=5,s=2,i=euclidean,n=1212){return super(t,s,i,n),super.parameter_list=["weight_adj","c"],this.parameter("weight_adj",e),this.parameter("c",r),this}
/**
     * 
     * @param {Matrix} [pca = null] - Initial Embedding (if null then PCA gets used). 
     * @param {KNN} [knn = null] - KNN Object (if null then BallTree gets used). 
     */init(t=null,e=null){const r=this.X,s=r.shape[0],i=this._d,n=this._metric,o=this._c;this.n_inliers=2*o,this.n_outliers=1*o,this.n_random=1*o,this.Y=t||new PCA(r,i).transform(),//.mult(.01);
this.knn=e||new BallTree(r.to2dArray,n);const{triplets:a,weights:h}=this._generate_triplets(this.n_inliers,this.n_outliers,this.n_random);return this.triplets=a,this.weights=h,this.lr=1e3*s/a.shape[0],this.C=1/0,this.tol=1e-7,this.vel=new Matrix(s,i,0),this.gain=new Matrix(s,i,1),this}
/**
     * Generates {@link n_inliers} x {@link n_outliers} x {@link n_random} triplets.
     * @param {Number} n_inliers 
     * @param {Number} n_outliers 
     * @param {Number} n_random 
     */_generate_triplets(t,e,r){const s=this._metric,i=this._weight_adj,n=this.X,o=n.shape[0],a=this.knn,h=Math.min(t+20,o),l=new Matrix(o,h),_=new Matrix(o,h);for(let t=0;t<o;++t)a.search(n.row(t),h+1).raw_data().filter((t=>0!=t.value)).sort(((t,e)=>t.value-e.value)).forEach(((e,r)=>{l.set_entry(t,r,e.element.index),_.set_entry(t,r,e.value)}));
// scale parameter
const c=new Float64Array(o);for(let t=0;t<o;++t)c[t]=Math.max((_.entry(t,3)+_.entry(t,4)+_.entry(t,5)+_.entry(t,6))/4,1e-10);const u=this._find_p(_,c,l);let d=this._sample_knn_triplets(u,l,t,e),f=d.shape[0];const p=new Float64Array(f);for(let t=0;t<f;++t){const e=d.entry(t,0),r=d.entry(t,2);p[t]=s(n.row(e),n.row(r))}let m=this._find_weights(d,u,l,p,c);if(r>0){const{random_triplets:t,random_weights:e}=this._sample_random_triplets(n,r,c);d=d.concat(t,"vertical"),m=Float64Array.from([...m,...e])}f=d.shape[0];let y=-1/0;for(let t=0;t<f;++t)isNaN(m[t])&&(m[t]=0),y<m[t]&&(y=m[t]);let w=-1/0;for(let t=0;t<f;++t)m[t]/=y,m[t]+=1e-4,m[t]=Math.log(1+i*m[t]),w<m[t]&&(w=m[t]);for(let t=0;t<f;++t)m[t]/=w;return{triplets:d,weights:m}}
/**
     * Calculates the similarity matrix P
     * @private
     * @param {Matrix} knn_distances - matrix of pairwise knn distances
     * @param {Float64Array} sig - scaling factor for the distances
     * @param {Matrix} nbrs - nearest neighbors
     * @returns {Matrix} pairwise similarity matrix
     */_find_p(t,e,r){const[s,i]=t.shape;return new Matrix(s,i,((s,i)=>Math.exp(-(t.entry(s,i)**2)/e[s]/e[r.entry(s,i)])))}
/**
     * Sample nearest neighbors triplets based on the similarity values given in P.
     * @private
     * @param {Matrix} P - Matrix of pairwise similarities between each point and its neighbors given in matrix nbrs.
     * @param {Matrix} nbrs - Nearest neighbors indices for each point. The similarity values are given in matrix {@link P}. Row i corresponds to the i-th point.
     * @param {Number} n_inliers - Number of inlier points.
     * @param {Number} n_outliers - Number of outlier points.
     * 
     */_sample_knn_triplets(t,e,r,s){const i=e.shape[0],n=new Matrix(i*r*s,3);for(let o=0;o<i;++o){let a=o*r*s;const h=this.__argsort(t.row(o).map((t=>-t)));for(let t=0;t<r;++t){let r=t*s;const l=e.entry(o,h[t]),_=this._rejection_sample(s,i,h.slice(0,t+1));for(let t=0;t<s;++t){const e=a+r+t,s=_[t];n.set_entry(e,0,o),n.set_entry(e,1,l),n.set_entry(e,2,s)}}}return n}
/**
     * Should do the same as np.argsort()
     * @private
     * @param {Array} A 
     */__argsort(t){return t.map(((t,e)=>({d:t,i:e}))).sort(((t,e)=>t.d-e.d)).map((t=>t.i))}
/**
     * Samples {@link n_samples} integers from a given interval [0, {@link max_int}] while rejection the values that are in the {@link rejects}.
     * @private
     * @param {*} n_samples 
     * @param {*} max_int 
     * @param {*} rejects 
     */_rejection_sample(t,e,r){const s=this._randomizer,i=linspace(0,e-1).filter((t=>r.indexOf(t)<0));return s.choice(i,Math.min(t,i.length-2))}
/**
     * Calculates the weights for the sampled nearest neighbors triplets
     * @private
     * @param {Matrix} triplets - Sampled Triplets.
     * @param {Matrix} P - Pairwise similarity matrix.
     * @param {Matrix} nbrs - nearest Neighbors
     * @param {Float64Array} outlier_distances - Matrix of pairwise outlier distances
     * @param {Float64Array} sig - scaling factor for the distances.
     */_find_weights(t,e,r,s,i){const n=t.shape[0],o=new Float64Array(n);for(let a=0;a<n;++a){const n=t.entry(a,0),h=r.row(n).indexOf(t.entry(a,1)),l=e.entry(n,h);let _=Math.exp(-(s[a]**2)/(i[n]*i[t.entry(a,2)]));_<1e-20&&(_=1e-20),o[a]=l/_}return o}
/**
     * Sample uniformly ranom triplets
     * @private
     * @param {Matrix} X - Data matrix.
     * @param {Number} n_random - Number of random triplets per point
     * @param {Float64Array} sig - Scaling factor for the distances
     */_sample_random_triplets(t,e,r){const s=this._metric,i=this._randomizer,n=t.shape[0],o=new Matrix(n*e,3),a=new Float64Array(n*e);for(let h=0;h<n;++h){const l=h*e,_=[...linspace(0,h-1),...linspace(h+1,n-1)];for(let n=0;n<e;++n){let[e,c]=i.choice(_,2),u=Math.exp(-(s(t.row(h),t.row(e))**2)/(r[h]*r[e]));u<1e-20&&(u=1e-20);let d=Math.exp(-(s(t.row(h),t.row(c))**2)/(r[h]*r[c]));d<1e-20&&(d=1e-20),u<d&&([e,c]=[c,e],[u,d]=[d,u]);const f=l+n;o.set_entry(f,0,h),o.set_entry(f,1,e),o.set_entry(f,2,c),a[f]=u/d}}return{random_triplets:o,random_weights:a}}
/**
     * Computes the gradient for updating the embedding.
     * @param {Matrix} Y - The embedding
     */_grad(t){const e=this.n_inliers,r=this.n_outliers,s=this.triplets,i=this.weights,[n,o]=t.shape,a=s.shape[0],h=new Matrix(n,o,0);let l=new Array(o).fill(0),_=new Array(o).fill(0),c=1,u=1,d=0,f=0;const p=n*e*r;for(let e=0;e<a;++e){const[n,a,m]=s.row(e);
// update y_ij, y_ik, d_ij, d_ik
if(e%r==0||e>=p){c=1,u=1;for(let e=0;e<o;++e){const r=t.entry(n,e),s=t.entry(a,e),i=t.entry(m,e);l[e]=r-s,_[e]=r-i,c+=l[e]**2,u+=_[e]**2}
// update y_ik and d_ik only
}else{u=1;for(let e=0;e<o;++e){const r=t.entry(n,e),s=t.entry(m,e);_[e]=r-s,u+=_[e]**2}}c>u&&++d,f+=i[e]/(1+u/c);const y=(i[e]/(c+u))**2;for(let t=0;t<o;++t){const e=l[t]*u*y,r=_[t]*c*y;h.set_entry(n,t,h.entry(n,t)+e-r),h.set_entry(a,t,h.entry(a,t)-e),h.set_entry(m,t,h.entry(m,t)+r)}}return{grad:h,loss:f,n_viol:d}}
/**
     * 
     * @param {Number} max_iteration 
     */transform(t=400){this.check_init();for(let e=0;e<t;++e)this._next(e);return this.projection}
/**
     * @yields {Matrix}
     * @returns {Matrix}
     */*generator(){this.check_init();for(let t=0;t<800;++t)this._next(t),yield this.projection;return this.projection}
/**
     * Does the iteration step.
     * @private
     * @param {Number} iter 
     */_next(t){const e=t>150?.5:.3,r=this.C,s=this.vel,i=this.Y.add(s.mult(e)),{grad:n,loss:o,n_viol:a}=this._grad(i);return this.C=o,this.Y=this._update_embedding(i,t,n),this.lr*=r>o+this.tol?1.01:.9,this.Y}
/**
     * Updates the embedding.
     * @private
     * @param {Matrix} Y 
     * @param {Number} iter 
     * @param {Matrix} grad 
     */_update_embedding(t,e,r){const[s,i]=t.shape,n=e>150?.9:.5,o=this.gain,a=this.vel,h=this.lr;for(let e=0;e<s;++e)for(let s=0;s<i;++s){const i=Math.sign(a.entry(e,s))!=Math.sign(r.entry(e,s))?o.entry(e,s)+.2:Math.max(.8*o.entry(e,s),.01);o.set_entry(e,s,i),a.set_entry(e,s,n*a.entry(e,s)-h*o.entry(e,s)*r.entry(e,s)),t.set_entry(e,s,t.entry(e,s)+a.entry(e,s))}return t}}
/**
 * @class
 * @alias Hierarchical_Clustering
 */class Hierarchical_Clustering{
/**
     * @constructor
     * @memberof module:clustering
     * @alias Hierarchical_Clustering
     * @todo needs restructuring. 
     * @param {Matrix} - Data or distance matrix if metric is 'precomputed'
     * @param {("single"|"complete"|"average")} [linkage = "complete"] 
     * @param {Function|"precomputed"} [metric = euclidean] 
     * @returns {Hierarchical_Clustering}
     */
constructor(t,e="complete",r=euclidean){if(this._id=0,this._matrix=t instanceof Matrix?t:Matrix.from(t),this._metric=r,this._linkage=e,"precomputed"===r&&t.shape[0]!==t.shape[1])throw"If metric is 'precomputed', then matrix has to be square!";return this.init(),this.root=this.do(),this}
/**
     * 
     * @param {Number} value - value where to cut the tree.
     * @param {("distance"|"depth")} [type = "distance"] - type of value.
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}.
     */get_clusters(t,e="distance"){let r,s=[];switch(e){case"distance":r=t=>t.dist;break;case"depth":r=t=>t.depth;break;default:throw"invalid type"}return this._traverse(this.root,r,t,s),s}
/**
     * @private
     * @param {} node 
     * @param {*} f 
     * @param {*} value 
     * @param {*} result 
     */_traverse(t,e,r,s){e(t)<=r?s.push(t.leaves()):(this._traverse(t.left,e,r,s),this._traverse(t.right,e,r,s))}
/**
     * computes the tree.
     */init(){const t=this._metric,e=this._matrix,r=this._n=e.shape[0],s=this._d_min=new Float64Array(r);let i;if("precomputed"!==t){i=new Matrix(r,r,0);//new Array(n);
for(let n=0;n<r;++n){s[n]=0;
//distance_matrix[i] = new Float64Array(n);
for(let o=0;o<r;++o)i.set_entry(n,o,n===o?1/0:t(e.row(n),e.row(o))),i.entry(n,s[n])>i.entry(n,o)&&(s[n]=o)}}else{i=this._matrix.clone();for(let t=0;t<r;++t)for(let e=0;e<r;++e)t===e?i.set_entry(t,e,1/0):i.entry(t,s[t])>i.entry(t,e)&&(s[t]=e)}this._distance_matrix=i;const n=this._clusters=new Array(r),o=this._c_size=new Uint16Array(r);for(let t=0;t<r;++t)n[t]=[],n[t][0]=new Cluster(this._id++,null,null,0,e.row(t),t,1,0),o[t]=1;return this}
/**
     * computes the tree.
     */do(){const t=this._n,e=this._d_min,r=this._distance_matrix,s=this._clusters,i=this._c_size,n=this._linkage;let o=null;for(let a=0,h=t-1;a<h;++a){let a=0;for(let s=0;s<t;++s)r.entry(s,e[s])<r.entry(a,e[a])&&(a=s);let h=e[a],l=s[a][0],_=s[h][0],c=l.isLeaf?[l.index]:l.index,u=_.isLeaf?[_.index]:_.index,d=c.concat(u),f=new Cluster(this._id++,l,_,r.entry(a,h),null,d);l.parent=f,_.parent=f,s[a].unshift(f),i[a]+=i[h];for(let e=0;e<t;++e)switch(n){case"single":r.entry(a,e)>r.entry(h,e)&&(r.set_entry(e,a,r.entry(h,e)),r.set_entry(a,e,r.entry(h,e)));break;case"complete":r.entry(a,e)<r.entry(h,e)&&(r.set_entry(e,a,r.entry(h,e)),r.set_entry(a,e,r.entry(h,e)));break;case"average":const t=(i[a]*r.entry(a,e)+i[h]*r.entry(h,e))/(i[a]+i[e]);r.set_entry(e,a,t),r.set_entry(a,e,t)}r.set_entry(a,a,1/0);for(let e=0;e<t;++e)r.set_entry(e,h,1/0),r.set_entry(h,e,1/0);for(let s=0;s<t;++s)e[s]===h&&(e[s]=a),r.entry(a,s)<r.entry(a,e[a])&&(e[a]=s);o=f}return o}}class Cluster{constructor(t,e,r,s,i,n,o,a){return this.id=t,this.left=e,this.right=r,this.dist=s,this.index=n,this.size=o??e.size+r.size,this.depth=a??1+Math.max(e.depth,r.depth),this.centroid=i??this._calculate_centroid(e,r),this.parent=null,this}_calculate_centroid(t,e){const r=t.size,s=e.size,i=t.centroid,n=e.centroid,o=this.size,a=t.centroid.length,h=new Float64Array(a);for(let t=0;t<a;++t)h[t]=(r*i[t]+s*n[t])/o;return h}get isLeaf(){return 0===this.depth}leaves(){if(this.isLeaf)return[this];const t=this.left,e=this.right;return(t.isLeaf?[t]:t.leaves()).concat(e.isLeaf?[e]:e.leaves())}descendants(){if(this.isLeaf)return[this];const t=this.left.descendants(),e=this.right.descendants();return t.concat(e).concat([this])}}
/**
 * @class
 * @alias KMeans
 */class KMeans{
/**
     * @constructor
     * @memberof module:clustering
     * @alias KMeans
     * @todo needs restructuring. 
     * @param {Matrix} matrix 
     * @param {Numbers} K 
     * @param {Function} [metric = euclidean] 
     * @param {Number} [seed = 1987]
     * @param {Boolean} [init = true]
     * @returns {KMeans}
     */
constructor(t,e,r=euclidean,s=1987,i=!0){this._metric=r,this._matrix=t,this._K=e;const[n,o]=t.shape;return this._N=n,this._D=o,e>n&&(e=n),this._randomizer=new Randomizer(s),this._clusters=new Array(n).fill(void 0),this._cluster_centroids=this._get_random_centroids(e),i&&this.init(e,this._cluster_centroids),this}
/**
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}. 
     */get_clusters(){const t=this._K,e=this._clusters,r=new Array(t).fill().map((()=>new Array));return e.forEach(((t,e)=>r[t].push(e))),r}
/**
     * @private
     * @param {Array} points 
     * @param {Array} candidates 
     */_furthest_point(t,e){const r=this._matrix,s=this._metric;let i=t.length;return Heap.heapify(e,(e=>{const n=r.row(e);let o=0;for(let e=0;e<i;++e)o+=s(n,t[e]);return o}),"max").pop().element}_get_random_centroids(t){const e=this._N,r=this._randomizer,s=this._matrix,i=new Array(t).fill(),n=linspace(0,e-1),o=r.random_int%(e-1);i[0]=s.row(o);const a=[o],h=Math.floor((e-t)/t);// / K
for(let e=1;e<t;++e){
// sampling + kmeans++ improvement?
const t=r.choice(n.filter((t=>-1==a.indexOf(t))),h),o=this._furthest_point(i.slice(0,e),t);a.push(o),i[e]=s.row(o)}return i}_iteration(t){const e=t.length,r=this._N,s=this._D,i=this._matrix,n=this._metric,o=this._clusters;let a=!1;
// find nearest cluster centroid.
for(let s=0;s<r;++s){const r=i.row(s);let h=1/0,l=null;for(let s=0;s<e;++s){let e=n(t[s],r);e<h&&(h=e,l=s)}o[s]!==l&&(a=!0),o[s]=l}
// update cluster centroid
// reset cluster centroids to 0
for(let r=0;r<e;++r){const e=t[r];for(let t=0;t<s;++t)e[t]=0}
// compute centroid
return this._compute_centroid(t),{clusters_changed:a,cluster_centroids:t}}_compute_centroid(t){const e=t.length,r=this._N,s=this._D,i=this._matrix,n=this._clusters,o=new Array(e).fill(0);for(let e=0;e<r;++e){const r=i.row(e),a=n[e];o[a]++;const h=t[a];for(let t=0;t<s;++t)h[t]+=r[t]}for(let r=0;r<e;++r){const e=o[r];t[r]=t[r].map((t=>t/e))}}
/**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */init(t,e){t||(t=this._K),e||(e=this._get_random_centroids(t));let r=!1;do{const t=this._iteration(e);e=t.cluster_centroids,r=t.clusters_changed}while(r)}}
/**
 * @class
 * @alias KMedoids
 */class KMedoids{
/**
     * @constructor
     * @memberof module:clustering
     * @alias KMedoids
     * @todo needs restructuring. 
     * @param {Matrix} matrix - data matrix
     * @param {Numbers} K - number of clusters
     * @param {number} [max_iter=null] - maximum number of iterations. Default is 10 * Math.log10(N)
     * @param {Function} [metric = euclidean] - metric defining the dissimilarity 
     * @param {Number} [seed = 1212] - seed value for random number generator
     * @returns {KMedoids}
     * @see {@link https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16} Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms
     */
constructor(t,e,r=null,s=euclidean,i=1212){this._metric=s,this._matrix=t,this._A=this._matrix.to2dArray,this._K=e;const[n,o]=t.shape;return this._N=n,this._D=o,this._max_iter=r||10*Math.log10(n),this._distance_matrix=new Matrix(n,n,"zeros"),
/* for (let i = 1; i < N; ++i) {
            for (let j = i + 1; j < N; ++j) {
                let dist = metric(this._A[i], this._A[j]);
                this._distance_matrix.set_entry(i, j, dist);
                this._distance_matrix.set_entry(j, i, dist)
            }
        } */
e>n&&(e=n),this._randomizer=new Randomizer(i),this._clusters=new Array(n).fill(void 0),this._cluster_medoids=this._get_random_medoids(e),
//if (init) this.init(K, this._cluster_medoids);
this._is_initialized=!1,this}
/**
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}. 
     */get_clusters(){const t=this._K,e=this._A;this._is_initialized||this.init(t,this._cluster_medoids);const r=new Array(t).fill().map((()=>new Array));return e.forEach(((t,e)=>{r[this._nearest_medoid(t,e).index_nearest].push(e)})),r.medoids=this._cluster_medoids,r}async*generator(){const t=this._max_iter;yield this.get_clusters();let e=!1,r=0;do{e=this._iteration(),yield this.get_clusters()}while(!e&&++r<t)}
/**
     * Algorithm 1. FastPAM1: Improved SWAP algorithm
     */
/* _iteration_1() {
        const A = this._A;
        const N = this._N;
        const K = this._K;
        const medoids = this._cluster_medoids;
        let DeltaTD = 0;
        let m0 = null;
        let x0 = null;
        A.forEach((x_j, j) => {
            if (medoids.findIndex(m => m === j) < 0) {
                const nearest_medoid = this._nearest_medoid(x_j, j);
                const d_j = nearest_medoid.distance_nearest; // distance to current medoid
                const deltaTD = new Array(K).fill(-d_j); // change if making j a medoid
                A.forEach((x_o, o) => {
                    // disance to new medoid
                    const d_oj = this._get_distance(o, j, x_o, x_j);
                    const {
                        "index_nearest": n,
                        "distance_nearest": d_n,
                        "distance_second": d_s,
                    } = this._nearest_medoid(x_o, o); 
                    this._clusters[o] = n; // cached values
                    deltaTD[n] += Math.min(d_oj, d_s) - d_n; // loss change
                    if (d_oj < d_n) { // reassignment check
                        deltaTD.forEach((d_i, i) => {
                            if (n !== i) {
                                deltaTD[i] = d_i + d_oj - d_n; // update loss change
                            }
                        });
                    }
                });
                // choose best medoid i;
                const i = deltaTD
                    .map((d, i) => [d, i])
                    .sort((d1, d2) => d1[0] - d2[0])[0][1];
                const deltaTD_i = deltaTD[i];
                // store
                if (deltaTD_i < DeltaTD) {
                    DeltaTD = deltaTD_i;
                    m0 = i;
                    x0 = j;
                }
            }
        });

        if (DeltaTD >= 0) {
            return true // break loop if DeltaTD >= 0
        }
        // swap roles of medoid m and non-medoid x;
        medoids[m0] = x0;
        this._cluster_medoids = medoids;
        return false
    } */
/** Algorithm 2. FastPAM2: SWAP with multiple candidates
     * 
     */_iteration(){const t=this._A,e=this._K,r=this._cluster_medoids,s=t.map(((t,e)=>this._nearest_medoid(t,e))),i=new Array(e).fill(0),n=new Array(e).fill(null);
// stop if no improvements were found
if(t.forEach(((o,a)=>{if(r.findIndex((t=>t===a))<0){const r=s[a].distance_nearest,h=new Array(e).fill(-r);// distance to current medoid
// change if making j a medoid
t.forEach(((t,r)=>{if(a===r)return;const i=this._get_distance(r,a,t,o),{index_nearest:n,distance_nearest:l,distance_second:_}=s[r];// distance to new medoid
// loss change for x_o
// Reassignment check
if(// cached
h[n]+=Math.min(i,_)-l,i<l)
// update loss change
for(let t=0;t<e;++t)t!==n&&(h[t]+=i-l)})),
// remember best swap for i;
h.map(((t,e)=>[t,e])).filter((([t,e])=>t<i[e])).forEach((([t,e])=>{t<i[e]&&(i[e]=t,n[e]=a)}))}})),Math.min(...i)>=0)return!0;
// execute all improvements
for(;Math.min(...i)<0;){
// swap roles of medoid m_i and non_medoid xs_i
const e=i.map(((t,e)=>[t,e])).sort((([t],[e])=>t-e))[0][1];0==r.filter((t=>t==n[e])).length&&(r[e]=n[e]),
// disable the swap just performed
i[e]=0,
// recompute TD for remaining swap candidates
i.map(((t,e)=>[t,e])).filter((([t])=>t<0)).forEach((([n,o])=>{const a=t[o];let h=0;t.forEach(((t,i)=>{r.findIndex((t=>t!=o&&t==i))>=0||e!=o&&(s[i].index_nearest===r[o]?h+=Math.min(this._get_distance(i,o,t,a),s[i].distance_second)-s[i].distance_nearest:h+=Math.min(this._get_distance(i,o,t,a)-s[i].distance_nearest,0))})),i[o]=h}))}return this._cluster_medoids=r,!1}_get_distance(t,e,r=null,s=null){if(t===e)return 0;const i=this._distance_matrix,n=this._A,o=this._metric;let a=i.entry(t,e);return 0===a&&(a=o(r||n[t],s||n[e]),i.set_entry(t,e,a),i.set_entry(e,t,a)),a}_nearest_medoid(t,e){const r=this._cluster_medoids,s=this._A,[i,n]=r.map(((r,i)=>{const n=s[r];return[this._get_distance(e,r,t,n),i]})).sort(((t,e)=>t[0]-e[0]));return{distance_nearest:i[0],index_nearest:i[1],distance_second:n[0],index_second:n[1]}}
/**
     * Computes {@link K} clusters out of the {@link matrix}.
     * @param {Number} K - number of clusters.
     */init(t,e){t||(t=this._K),e||(e=this._get_random_medoids(t));const r=this._max_iter;let s=!1,i=0;do{s=this._iteration()}while(!s&&++i<r);return this}
/**
     * Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.
     * @param {number} K - number of clusters
     * 
     */_get_random_medoids(t){const e=this._N,r=this._A,s=linspace(0,e-1),i=this._randomizer,n=Math.min(e,10+Math.ceil(Math.sqrt(e))),o=new Array(n).fill(1/0),a=[];
// first medoid
let h=1/0,l=i.choice(s,n);for(let t=0;t<n;++t){const e=l[t],s=r[e];for(let e=0;e<n;++e){if(e===t)continue;const i=r[l[e]];o[t]+=this._get_distance(t,e,s,i)}o[t]<h&&(h=o[t],// smallest distance sum
a.push(e))}
// other medoids
for(let e=1;e<t;++e){let t=1/0;l=i.choice(s.filter((t=>a.findIndex((e=>e===t))<0)),n);for(let e=0;e<n;++e){let s=0;const i=l[e],o=r[i];for(let t=0;t<n;++t){if(t===e)continue;const n=l[t],h=r[n];let _=this._get_distance(i,n,o,h)-Math.min(...a.map((t=>this._get_distance(n,t,h))));_<0&&(s+=_)}
// best reduction
s<t&&(t=s,a.push(i))}h+=t}return a.slice(0,t)}}
/**
 * @class
 * @alias OPTICS
 */class OPTICS{
/**
     * **O**rdering **P**oints **T**o **I**dentify the **C**lustering **S**tructure.
     * @constructor
     * @memberof module:clustering
     * @alias OPTICS
     * @todo needs restructuring. 
     * @param {Matrix} matrix - the data.
     * @param {Number} epsilon - the minimum distance which defines whether a point is a neighbor or not.
     * @param {Number} min_points - the minimum number of points which a point needs to create a cluster. (Should be higher than 1, else each point creates a cluster.)
     * @param {Function} [metric = euclidean] - the distance metric which defines the distance between two points of the {@link matrix}.
     * @returns {OPTICS}
     * @see {@link https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf}
     * @see {@link https://en.wikipedia.org/wiki/OPTICS_algorithm}
     */
constructor(t,e,r,s=euclidean){return this._matrix=t,this._epsilon=e,this._min_points=r,this._metric=s,this._ordered_list=[],this._clusters=[],this._DB=new Array(t.shape[0]).fill(),this.init(),this}
/**
     * Computes the clustering.
     */init(){const t=this._ordered_list,e=this._matrix,r=e.shape[0],s=this._DB,i=this._clusters;let n=this._cluster_index=0;for(let t=0;t<r;++t)s[t]={element:e.row(t),index:t,reachability_distance:void 0,processed:!1};for(const e of s)if(!e.processed&&(e.neighbors=this._get_neighbors(e),e.processed=!0,i.push([e.index]),n=i.length-1,t.push(e),null!=this._core_distance(e))){const t=new Heap(null,(t=>t.reachability_distance),"min");this._update(e,t),this._expand_cluster(t,i[n])}return this}
/**
     * 
     * @private
     * @param {Object} p - a point of {@link matrix}.
     * @returns {Array} An array consisting of the {@link epsilon}-neighborhood of {@link p}.
     */_get_neighbors(t){if("neighbors"in t)return t.neighbors;const e=this._DB,r=this._metric,s=this._epsilon,i=[];for(const n of e)n.index!=t.index&&r(t.element,n.element)<s&&i.push(n);return i}
/**
     * 
     * @private
     * @param {Object} p - a point of {@link matrix}.
     * @returns {Number} The distance to the {@link min_points}-th nearest point of {@link p}, or undefined if the {@link epsilon}-neighborhood has fewer elements than {@link min_points}.
     */_core_distance(t){const e=this._min_points,r=this._metric;if(!(t.neighbors&&t.neighbors.length<=e))return r(t.element,t.neighbors[e].element)}
/**
     * Updates the reachability distance of the points.
     * @private
     * @param {Object} p 
     * @param {Heap} seeds 
     */_update(t,e){const r=this._metric,s=this._core_distance(t),i=this._get_neighbors(t);//p.neighbors;
for(const n of i){if(n.processed)continue;const i=Math.max(s,r(t.element,n.element));
//if (q.reachability_distance == undefined) { // q is not in seeds
e.raw_data().findIndex((t=>t.element==n))<0?(n.reachability_distance=i,e.push(n)):// q is in seeds
i<n.reachability_distance&&(n.reachability_distance=i,e=Heap.heapify(e.data(),(t=>t.reachability_distance),"min"))}}
/**
     * Expands the {@link cluster} with points in {@link seeds}.
     * @private
     * @param {Heap} seeds 
     * @param {Array} cluster 
     */_expand_cluster(t,e){const r=this._ordered_list;for(;!t.empty;){const s=t.pop().element;s.neighbors=this._get_neighbors(s),s.processed=!0,e.push(s.index),r.push(s),null!=this._core_distance(s)&&(this._update(s,t),this._expand_cluster(t,e))}}
/**
     * Returns an array of clusters.
     * @returns {Array<Array>} Array of clusters with the indices of the rows in given {@link matrix}.
     */get_clusters(){const t=[],e=[],r=this._min_points;for(const s of this._clusters)s.length<r?e.push(...s):t.push(s);return t.push(e),t}
/**
     * @returns {Array} Returns an array, where the ith entry defines the cluster affirmation of the ith point of {@link matrix}. (-1 stands for outlier)
     */get_cluster_affirmation(){const t=this._matrix.shape[0],e=new Array(t).fill(),r=this.get_clusters();for(let t=0,s=r.length;t<s;++t){const i=r[t];for(const r of i)e[r]=t<s-1?t:-1}return e}}
/**
 * @class
 * @alias LSP
 */class LSP extends DR$1{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LSP
     * @param {Matrix} X - the high-dimensional data. 
     * @param {number} [k = Math.max(Math.floor(N / 10), 2)] - number of neighbors to consider.
     * @param {number} [control_points = Math.ceil(Math.sqrt(N))] - number of controlpoints
     * @param {number} [d = 2] - the dimensionality of the projection.
     * @param {function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @returns {LSP}
     * @see {@link https://ieeexplore.ieee.org/document/4378370}
     */
constructor(t,e,r,s=2,i=euclidean,n=1212){return super(t,s,i,n),super.parameter_list=["k","control_points"],this.parameter("k",Math.min(e??Math.max(Math.floor(this._N/10),2),this._N-1)),this.parameter("control_points",Math.min(r??Math.ceil(Math.sqrt(this._N)),this._N-1)),this._is_initialized=!1,this}
/**
     * 
     * @param {DR} DR - method used for position control points.
     * @param {DR_parameters} DR_parameters - array containing parameters for the DR method which projects the control points
     * @returns {LSP} 
     */init(t=MDS,e=[],r=BallTree){if(this._is_initialized)return this;const s=this.X,i=this._N,n=this.parameter("k"),o=this._d,a=this._metric,h=this.parameter("control_points"),l=new KMedoids(s,h,null,a).get_clusters().medoids,_=new Matrix(h,i,"zeros");l.forEach(((t,e)=>{_.set_entry(e,t,1)}));const c=new t(Matrix.from(l.map((t=>s.row(t)))),...e,o).transform(),u=s.to2dArray,d=new r(u,a),f=new Matrix(i,i,"I"),p=-1/n;u.forEach(((t,e)=>{for(const{index:r}of d.search(t,n).iterate())e!==r&&f.set_entry(e,r,p)}));const m=f.concat(_,"vertical"),y=new Matrix(i,o,"zeros").concat(c,"vertical");return this._A=m,this._b=y,this._is_initialized=!0,this}
/**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */transform(){this.check_init();const t=this._A,e=t.T,r=this._b,s=e.dot(t),i=e.dot(r);return this.Y=Matrix.solve_CG(s,i,this._randomizer),this.projection}}class TopoMap extends DR$1{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias Topomap
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {TopoMap}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
constructor(t,e=2,r=euclidean,s=1212){return super(t,e,r,s),super.parameter_list=[],[this._N,this._D]=this.X.shape,this._distance_matrix=new Matrix(this._N,this._N,0),this}
/**
     * @private
     */__lazy_distance_matrix(t,e,r){const s=this._distance_matrix,i=this.X,n=s.entry(t,e);if(0===n){let n=r(i.row(t),i.row(e));return s.set_entry(t,e,n),s.set_entry(e,t,n),n}return n}
/**
     * Computes the minimum spanning tree, using a given metric
     * @private
     * @param {Function} metric 
     * @see {@link https://en.wikipedia.org/wiki/Kruskal%27s_algorithm}
     */_make_minimum_spanning_tree(t=euclidean){const e=this._N,r=[...this.X];let s=new DisjointSet(r);const i=[];let n=[];for(let r=0;r<e;++r)for(let s=r+1;s<e;++s)n.push([r,s,this.__lazy_distance_matrix(r,s,t)]);n=n.sort(((t,e)=>t[2]-e[2]));for(const[t,e,o]of n){const n=s.find(r[t]),a=s.find(r[e]);n!==a&&(i.push([t,e,o]),s.union(n,a))}return i.sort(((t,e)=>t[2]-e[2]))}
/**
     * initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree.
     */init(){return this.Y=new Matrix(this._N,this._d,0),this._Emst=this._make_minimum_spanning_tree(this._metric),this._is_initialized=!0,this}
/**
     * Returns true if Point C is left of line AB.
     * @private
     * @param {Array} PointA - Point A of line AB
     * @param {Array} PointB - Point B of line AB
     * @param {Array} PointC - Point C
     * @returns {Boolean}
     */__hull_cross([t,e],[r,s],[i,n]){return(r-t)*(n-e)-(s-e)*(i-t)<=0}
/**
     * Computes the convex hull of the set of Points S
     * @private
     * @param {Array} S - Set of Points.
     * @see {@link https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#JavaScript}
     * @returns {Array} convex hull of S. Starts at the bottom-most point and continues counter-clockwise.
     */__hull(t){const e=t.sort((([t,e],[r,s])=>e-s||t-r)),r=e.length;if(r<=2)return e;const s=[];for(let t=0;t<r;++t){for(;s.length>=2&&this.__hull_cross(s[s.length-2],s[s.length-1],e[t]);)s.pop();s.push(e[t])}const i=[];for(let t=r-1;t>=0;--t){for(;i.length>=2&&this.__hull_cross(i[i.length-2],i[i.length-1],e[t]);)i.pop();i.push(e[t])}return i.pop(),s.pop(),s.concat(i)}
/**
     * Finds the angle to rotate Point A and B to lie on a line parallel to the x-axis.
     * @private
     * @param {Array} PointA 
     * @param {Array} PointB
     * @return {Object} Object containing the sinus- and cosinus-values for a rotation.  
     */__findAngle([t,e],[r,s]){const i=euclidean([t,e],[r,s]);if(0===i)return{sin:0,cos:1};const n=[(r-t)/i,(s-e)/i],o=n[0];let a=Math.sqrt(1-o*o);return a=n[1]>=0?-a:a,{sin:a,cos:o}}
/**
     * @private
     * @param {Array} hull 
     * @param {Array} p 
     * @param {Bool} topEdge
     */__align_hull(t,e,r){let s,i,n,o=-1;for(let r=0;r<t.length;++r){const i=euclidean(t[r],e);(-1===o||s>i)&&(s=i,o=r)}r?(i=t[o],n=t[(o+1)%t.length]):(0==o&&(o=t.length-1),i=t[o],n=t[(o-1)%t.length]);const a={tx:-t[o][0],ty:-t[o][1]};if(t.length>=2){const{sin:t,cos:e}=this.__findAngle(i,n);a.sin=t,a.cos=e}else a.sin=0,a.cos=1;return a}
/**
     * @private
     * @param {Array} Point - The point which should get transformed.
     * @param {Object} Transformation - contains the values for translation and rotation.
     */__transform([t,e],{tx:r,ty:s,sin:i,cos:n}){let o=t+r,a=e+s;return[o*n-a*i,o*i+a*n]}
/**
     * Calls {@link __transform} for each point in Set C
     * @private
     * @param {Array} C - Set of points.
     * @param {Object} t - Transform object. 
     * @param {Number} yOffset - value to offset set C.
     */__transform_component(t,e,r){const s=t.length;for(let i=0;i<s;++i){const s=t[i],[n,o]=this.__transform(s,e);s[0]=n,s[1]=o+r}}
/**
     * @private
     * @param {Array} u - point u
     * @param {Array} v - point v
     * @param {Number} w - edge weight w
     */__align_components(t,e,r){const s=[...t.__disjoint_set.children],i=[...e.__disjoint_set.children],n=this.__hull(s),o=this.__hull(i),a=this.__align_hull(n,t,!1),h=this.__align_hull(o,e,!0);this.__transform_component(s,a,0),this.__transform_component(i,h,r)}
/**
     * Transforms the inputdata {@link X} to dimensionality 2.
     */transform(){this._is_initialized||this.init();const t=this._Emst,e=[...this.Y],r=new DisjointSet(e.map(((t,e)=>(t.i=e,t))));for(const[s,i,n]of t){const t=r.find(e[s]),o=r.find(e[i]);t!==o&&(this.__align_components(t,o,n),r.union(t,o))}return this.projection}*generator(){this._is_initialized||this.init();const t=this._Emst,e=[...this.Y],r=new DisjointSet(e.map(((t,e)=>(t.i=e,t))));for(const[s,i,n]of t){const t=r.find(e[s]),o=r.find(e[i]);t!==o&&(this.__align_components(t,o,n),r.union(t,o),yield this.projection)}return this.projection}}
/**
 * @class
 * @alias SAMMON
 */class SAMMON extends DR{
/**
     * 
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias SAMMON
     * @param {Matrix} X - the high-dimensional data. 
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.  
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @returns {SAMMON}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
constructor(t,e=.1,r=2,s=euclidean,i=1212){return super(t,r,s,i),super.parameter_list=["magic"],this.parameter("magic",e),[this._N,this._D]=this.X.shape,this}
/**
     * initializes SAMMON. Sets all projcted points to zero, and computes a minimum spanning tree.
     */init(t="random",e=null){const r=this._N,s=this._d;if("random"===t){const t=this._randomizer;this.Y=new Matrix(r,s,(()=>t.random))}else t instanceof DR$1&&(this.Y=t.transform(this.X));return this.distance_matrix=e||this.__distance_matrix(this.X),this}
/**
     * @private
     * @param {Matrix} A
     * @returns {Matrix} 
     */__distance_matrix(t){const e=this._metric,r=t.shape[0],s=new Matrix(r,r);for(let i=0;i<r;++i){const n=t.row(i);for(let o=i;o<r;++o){let r=i===o?0:e(n,t.row(o));s.set_entry(i,o,r),s.set_entry(o,i,r)}}return s}
/**
     * Transforms the inputdata {@link X} to dimenionality 2.
     */transform(t=200){this._is_initialized||this.init();for(let e=0;e<t;++e)this._step();return this.projection}*generator(t=200){this._is_initialized||this.init();for(let e=0;e<t;++e)this._step(),yield this.projection;return this.projection}_step(){const t=this.parameter("magic"),e=this.distance_matrix,r=this._N,s=this._d,i=this._metric;let n=this.Y,o=new Matrix(r,s,0),a=new Float64Array(s);for(let h=0;h<r;++h){let l=new Float64Array(s),_=new Float64Array(s);const c=n.row(h);for(let t=0;t<r;++t){if(h===t)continue;const r=n.row(t),o=new Float64Array(s);for(let t=0;t<s;++t)o[t]=c[t]-r[t];const a=i(c,r),u=e.entry(h,t),d=u-a,f=Math.max(u*a,.01);for(let t=0;t<s;++t)l[t]+=o[t]*d/f,_[t]+=(d-Math.pow(o[t],2)*(1+d/a)/a)/f}for(let e=0;e<s;++e){const r=n.entry(h,e)+(t*l[e]/Math.abs(_[e])||0);o.set_entry(h,e,r),a[e]+=r}}for(let t=0;t<s;++t)a[t]/=r;for(let t=0;t<r;++t)for(let e=0;e<s;++e)n.set_entry(t,e,o.entry(t,e)-a[e]);return n}}var t="0.3.16";export{BallTree,DisjointSet,FASTMAP,Heap,Hierarchical_Clustering,ISOMAP,KMeans,KMedoids,KNN,LDA,LLE,LSP,LTSA,MDS,Matrix,OPTICS,PCA,Randomizer,SAMMON,TSNE,TopoMap,TriMap,UMAP,canberra,chebyshev,cosine,distance_matrix,euclidean,euclidean_squared,hamming,jaccard,k_nearest_neighbors,kahan_sum,linspace,manhattan,max,neumair_sum,norm,powell,qr,simultaneous_poweriteration$1 as simultaneous_poweriteration,sokal_michener,t as version,yule};
