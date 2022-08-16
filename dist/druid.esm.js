// https://renecutura.eu v0.6.3 Copyright 2022 Rene Cutura
/**
 * Computes the euclidean distance (<code>l<sub>2</sub></code>) between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias euclidean
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the euclidean distance between <code>a</code> and <code>b</code>.
 */
function euclidean(t,e){return Math.sqrt(euclidean_squared(t,e))}
/**
 * Computes the squared euclidean distance (l<sub>2</sub><sup>2</sup>) between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias euclidean_squared
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the squared euclidean distance between <code>a</code> and <code>b</code>.
 */function euclidean_squared(t,e){if(t.length!=e.length)return;const r=t.length;let s=0;for(let i=0;i<r;++i){const r=t[i]-e[i];s+=r*r}return s}
/**
 * Computes the cosine distance (not similarity) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias cosine
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} The cosine distance between {@link a} and {@link b}.
 * 
 * @example
 * import * as druid from "@saehrimnir/druidjs";
 * 
 * druid.cosine([1,0],[1,1]) == 0.7853981633974484 == π/4;
 * 
 */function cosine(t,e){if(t.length!==e.length)return;let r=t.length,s=0,i=0,n=0;for(let o=0;o<r;++o)s+=t[o]*e[o],i+=t[o]*t[o],n+=e[o]*e[o];return Math.acos(s/(Math.sqrt(i)*Math.sqrt(n)))}
/**
 * Computes the manhattan distance (<code>l<sub>1</sub></code>) between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias manhattan
 * @param {Array<Number>} a
 * @param {Array<Number>} b
 * @returns {Number} the manhattan distance between <code>a</code> and <code>b</code>.
 */function manhattan(t,e){if(t.length!=e.length)return;const r=t.length;let s=0;for(let i=0;i<r;++i)s+=Math.abs(t[i]-e[i]);return s}
/**
 * Computes the chebyshev distance (L<sub>∞</sub>) between {@link a} and {@link b}.
 * @memberof module:metrics
 * @alias chebyshev
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the chebyshev distance between {@link a} and {@link b}.
 */function chebyshev(t,e){if(t.length!=e.length)return;const r=t.length;let s=[];for(let i=0;i<r;++i)s.push(Math.abs(t[i]-e[i]));return Math.max(...s)}
/**
 * Computes the canberra distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias canberra
 * @param {Number[]} a 
 * @param {Number[]} b 
 * @returns {Number} the canberra distance between <code>a</code> and <code>b</code>.
 * @see {@link https://en.wikipedia.org/wiki/Canberra_distance}
 */function canberra(t,e){if(t.length!==e.length)return;const r=t.length;let s=0;for(let i=0;i<r;++i)s+=Math.abs(t[i]-e[i])/(Math.abs(t[i])+Math.abs(e[i]));return s}
/**
 * Computes the jaccard distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias jaccard
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the jaccard distance between <code>a</code> and <code>b</code>.
 */function jaccard(t,e){if(t.length!=e.length)return;const r=t.length;let s=0,i=0;for(let n=0;n<r;++n){const r=0!=t[n],o=0!=e[n];s+=r||o,i+=r&&o}return(s-i)/s}
/**
 * Computes the hamming distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias hamming
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the hamming distance between <code>a</code> and <code>b</code>.
 */function hamming(t,e){if(t.length!=e.length)return;const r=t.length;let s=0;for(let i=0;i<r;++i){s+=t[i]!=e[i]}return s/r}
/**
 * Computes the Sokal-Michener distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias sokal_michener
 * @param {Number[]} a 
 * @param {Number[]} b 
 * @returns {Number} the Sokal-Michener distance between <code>a</code> and <code>b</code>.  
 */function sokal_michener(t,e){if(t.length!=e.length)return;const r=t.length;let s=0;for(let i=0;i<r;++i){s+=0!=t[i]!=(0!=e[i])}return 2*s/(r+s)}
/**
 * Computes the yule distance between <code>a</code> and <code>b</code>.
 * @memberof module:metrics
 * @alias yule
 * @param {Number[]} a
 * @param {Number[]} b
 * @returns {Number} the yule distance between <code>a</code> and <code>b</code>.
 */function yule(t,e){if(t.length!=e.length)return;const r=t.length;let s=0,i=0,n=0;for(let o=0;o<r;++o){const r=0!=t[o],a=0!=e[o];s+=r&&a,i+=r&&!a,n+=!r&&r}return 0==i||0==n?0:2*i*n/(s*(r-s-i-n)+i*n)}
/**
 * Computes the k-nearest neighbors of each row of {@link A}.
 * @memberof module:matrix
 * @alias k_nearest_neigbhors
 * @param {Matrix} A - Either the data matrix, or a distance matrix.
 * @param {Number} k - The number of neighbors to compute.
 * @param {Function|"precomputed"} [metric=euclidean]
 * @returns {Array<Object>} -
 */function k_nearest_neighbors(t,e,r=euclidean){const s=t.shape[0];let i="precomputed"==r?t:distance_matrix(t,r),n=new Array(s);for(let t=0;t<s;++t)n[t]=Array.from(i.row(t)).map(((e,r)=>({i:t,j:r,distance:e}))).sort(((t,e)=>t.distance-e.distance)).slice(1,e+1);return n}
/**
 * Computes the distance matrix of datamatrix {@link A}.
 * @memberof module:matrix
 * @alias distance_matrix
 * @param {Matrix} A - Matrix.
 * @param {Function} [metric=euclidean] - The diistance metric.
 * @returns {Matrix} D - The distance matrix of {@link A}.
 */function distance_matrix(t,e=euclidean){let r=t.shape[0];const s=new Matrix(r,r);for(let i=0;i<r;++i){const n=t.row(i);for(let o=i+1;o<r;++o){const r=e(n,t.row(o));s.set_entry(i,o,r),s.set_entry(o,i,r)}}return s}
/**
 * Creates an Array containing {@link number} numbers from {@link start} to {@link end}.
 * If <code>{@link number} = null</null>.
 * @memberof module:matrix
 * @alias linspace
 * @param {Number} start - Start value.
 * @param {Number} end - End value.
 * @param {Number} [number = null] - Number of number between {@link start} and {@link end}.
 * @returns {Array} - An array with {@link number} entries, beginning at {@link start} ending at {@link end}.
 */function linspace(t,e,r=null){if(r||(r=Math.max(Math.round(e-t)+1,1)),r<2)return 1===r?[t]:[];let s=new Array(r);for(let i=r-=1;i>=0;--i)s[i]=(i*e+(r-i)*t)/r;return s}
//import { neumair_sum } from "../numerical/index";
/**
 * Computes the norm of a vector, by computing its distance to **0**.
 * @memberof module:matrix
 * @alias norm
 * @param {Matrix|Array<Number>|Float64Array} v - Vector.
 * @param {Function} [metric = euclidean] - Which metric should be used to compute the norm.
 * @returns {Number} - The norm of {@link v}.
 */function norm(t,e=euclidean){let r=null;if(t instanceof Matrix){let[e,s]=t.shape;if(1===e)r=t.row(0);else{if(1!==s)throw new Error("Matrix must be 1d!");r=t.col(0)}}else r=t;const s=r.length;return e(r,new Float64Array(s))}
/**
 * Normalizes Vector {@link v}.
 * @memberof module:matrix
 * @alias normalize
 * @param {Array<Number>|Float64Array} v - Vector
 * @param {Function} metric 
 * @returns {Array<Number>|Float64Array} - The normalized vector with length 1.
 */function normalize(t,e=euclidean){const r=norm(t,e);return t.map((t=>t/r))}
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
 * @param {Number[]} summands - Array of values to sum up.
 * @returns {Number} The sum.
 * @see {@link https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements}
 */function neumair_sum(t){const e=t.length;let r=0,s=0;for(let i=0;i<e;++i){const e=t[i],n=r+e;Math.abs(r)>=Math.abs(e)?s+=r-n+e:s+=e-n+r,r=n}return r+s}
/**
 * Computes the QR Decomposition of the Matrix `A` using Gram-Schmidt process.
 * @memberof module:linear_algebra
 * @alias qr
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_the_Gram%E2%80%93Schmidt_process}
 */function qr_gramschmidt(t){const[e,r]=t.shape,s=new Matrix(e,r,"identity"),i=new Matrix(r,r,0);for(let n=0;n<r;++n){let r=t.col(n);for(let t=0;t<n;++t){const o=s.col(t),a=neumair_sum(o.map(((t,e)=>t*r[e])));for(let t=0;t<e;++t)r[t]-=a*o[t];i.set_entry(t,n,a)}const o=norm(r,euclidean);for(let t=0;t<e;++t)s.set_entry(t,n,r[t]/o);i.set_entry(n,n,o)}return{R:i,Q:s}}
/**
 * Computes the QR Decomposition of the Matrix {@link A} with householder transformations.
 * @memberof module:linear_algebra
 * @alias qr_householder
 * @param {Matrix} A
 * @returns {{R: Matrix, Q: Matrix}}
 * @see {@link https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections}
 * @see {@link http://mlwiki.org/index.php/Householder_Transformation}
 */function qr_householder(t){const[e,r]=t.shape,s=new Matrix(e,e,"I"),i=t.clone();for(let t=0;t<r;++t){const e=Matrix.from(i.col(t).slice(t)),r=norm(e),n=e.entry(0,0),o=-Math.sign(n),a=n-o*r,h=e.divide(a).set_entry(0,0,1),l=-o*a/r,_=h.outer(h),c=i.get_block(t,0),u=c.sub(_.dot(c).mult(l)),d=s.get_block(0,t),m=d.sub(d.dot(_).mult(l));i.set_block(t,0,u),s.set_block(0,t,m)}return{R:i,Q:s}}
/**
 * Computes the `k` biggest Eigenvectors and Eigenvalues from Matrix `A` with the QR-Algorithm.
 * @memberof module:linear_algebra
 * @alias simultaneous_poweriteration
 * @param {Matrix} A - The Matrix
 * @param {Number} k - The number of eigenvectors and eigenvalues to compute.
 * @param {Object} parameters - Object containing parameterization of the simultanious poweriteration method.
 * @param {Number} [parameters.max_iterations=100] - The number of maxiumum iterations the algorithm should run.
 * @param {Number|Randomizer} [parameters.seed=1212] - The seed value or a randomizer used in the algorithm.
 * @param {Function} [parameters.qr=qr_gramschmidt] - The QR technique to use.
 * @param {Number} [parameters.tol=1e-8] - Tolerated error for stopping criteria.
 * @returns {{eigenvalues: Number[], eigenvectors: Number[][]}} the `k` biggest eigenvectors and eigenvalues of Matrix `A`.
 */function simultaneous_poweriteration(t,e=2,{seed:r=1212,max_iterations:s=100,qr:i=qr_gramschmidt,tol:n=1e-8}={}){const o=r instanceof Randomizer?r:new Randomizer(r);t instanceof Matrix||(t=Matrix.from(t));const a=t.shape[0];let{Q:h,R:l}=i(new Matrix(a,e,(()=>2*(o.random-.5))));for(;s--;){const e=h,r=i(t.dot(h));h=r.Q,l=r.R;if(euclidean_squared(h.values,e.values)<n)break}return{eigenvalues:l.diag,eigenvectors:h.transpose().to2dArray}}
/**
 * Computes the inner product between two arrays of the same length.
 * @memberof module:linear_algebra
 * @alias inner_product
 * @param {Array|Float64Array} a - Array a
 * @param {Array|Float64Array} b - Array b
 * @returns The inner product between {@link a} and {@link b}
 */function inner_product(t,e){const r=t.length;if(r!=e.length)throw new Error("Array a and b must have the same length!");let s=0;for(let i=0;i<r;++i)s+=t*e;return s}
/**
 * @class
 * @alias Matrix
 * @requires module:numerical/neumair_sum
 */class Matrix{
/**
     * creates a new Matrix. Entries are stored in a Float64Array.
     * @memberof module:matrix
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
     * let S = Matrix.from([1, 2, 3], "diag"); // creates a 3 by 3 matrix with 1, 2, 3 on its diagonal. [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
     */static from(t,e="row"){if(t instanceof Matrix)return t.clone();if(Matrix.isArray(t)){let r=t.length;if(0===r)throw new Error("Array is empty");
// 1d
if(Matrix.isArray(t[0])){let e=t[0].length;for(let s=0;s<r;++s)if(t[s].length!==e)throw new Error("various array lengths");return new Matrix(r,e,((e,r)=>t[e][r]))}if("row"===e)return new Matrix(1,r,((e,r)=>t[r]));
// 2d
if("col"===e)return new Matrix(r,1,(e=>t[e]));if("diag"===e)return new Matrix(r,r,((e,r)=>e==r?t[e]:0));throw new Error("1d array has NaN entries")}if("number"==typeof t)return new Matrix(1,1,t);throw new Error("error")}
/**
     * Returns the {@link row}<sup>th</sup> row from the Matrix.
     * @param {Number} row
     * @returns {Float64Array}
     */row(t){const e=this.values,r=this._cols;return e.subarray(t*r,(t+1)*r)}
/**
     * Returns an generator yielding each row of the Matrix.
     * @yields {Float64Array}
     */*iterate_rows(){const t=this._cols,e=this._rows,r=this.values;for(let s=0;s<e;++s)yield r.subarray(s*t,(s+1)*t)}
/**
     * Makes a {@link Matrix} object an iterable object.
     * @yields {Float64Array}
     */*[Symbol.iterator](){for(const t of this.iterate_rows())yield t}
/**
     * Sets the entries of {@link row}<sup>th</sup> row from the Matrix to the entries from {@link values}.
     * @param {Number} row
     * @param {Array} values
     * @returns {Matrix}
     */set_row(t,e){const r=this._cols;if(Matrix.isArray(e)&&e.length===r){const s=t*r;for(let t=0;t<r;++t)this.values[s+t]=e[t]}else{if(!(e instanceof Matrix&&e.shape[1]===r&&1===e.shape[0]))throw new Error("Values not valid! Needs to be either an Array, a Float64Array, or a fitting Matrix!");{const s=t*r;for(let t=0;t<r;++t)this.values[s+t]=e._data[t]}}return this}
/**
     * Swaps the rows {@link row1} and {@link row2} of the Matrix.
     * @param {Number} row1
     * @param {Number} row2
     * @returns {Matrix}
     */swap_rows(t,e){const r=this._cols,s=this.values;for(let i=t*r,n=e*r,o=0;o<r;++o,++i,++n){const t=s[i];s[i]=s[n],s[n]=t}}
/**
     * Returns the {@link col}<sup>th</sup> column from the Matrix.
     * @param {Number} col
     * @returns {Array}
     */col(t){const e=new Float64Array(this._rows);for(let r=0;r<this._rows;++r)e[r]=this.values[r*this._cols+t];return e}
/**
     * Returns the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix.
     * @param {int} row
     * @param {int} col
     * @returns {float64}
     */entry(t,e){return this.values[t*this._cols+e]}
/**
     * Sets the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix to the given {@link value}.
     * @param {int} row
     * @param {int} col
     * @param {float64} value
     * @returns {Matrix}
     */set_entry(t,e,r){return this.values[t*this._cols+e]=r,this}
/**
     * Adds a given {@link value} to the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix.
     * @param {int} row
     * @param {int} col
     * @param {float64} value
     * @returns {Matrix}
     */add_entry(t,e,r){return this.values[t*this._cols+e]+=r,this}
/**
     * Subtracts a given {@link value} from the {@link col}<sup>th</sup> entry from the {@link row}<sup>th</sup> row of the Matrix.
     * @param {int} row
     * @param {int} col
     * @param {float64} value
     * @returns {Matrix}
     */sub_entry(t,e,r){return this.values[t*this._cols+e]-=r,this}
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
     */inverse(){const t=this._rows,e=this._cols,r=this.clone(),s=new Matrix(t,e,"I");
// foreach column
for(let i=0;i<e;++i){
// Search for maximum in this column (pivot)
let n=i,o=Math.abs(r.entry(i,i));for(let e=i+1;e<t;++e){const t=Math.abs(r.entry(e,i));o<t&&(n=e,o=t)}if(0===o)throw new Error("Cannot compute inverse of Matrix, determinant is zero");
// Swap maximum row with current row
n!==i&&(r.swap_rows(i,n),s.swap_rows(i,n));
// eliminate non-zero values on the other rows at column c
const a=r.row(i),h=s.row(i);for(let n=0;n<t;++n)if(n!==i){
// eliminate value at column c and row r
const t=r.row(n),o=s.row(n);if(0!==t[i]){const r=t[i]/a[i];
// sub (f * row c) from row r to eliminate the value at column c
for(let s=i;s<e;++s)t[s]-=r*a[s];for(let t=0;t<e;++t)o[t]-=r*h[t]}}else{
// normalize value at Acc to 1 (diagonal):
// divide each value of row r=c by the value at Acc
const t=a[i];for(let r=i;r<e;++r)a[r]/=t;for(let r=0;r<e;++r)h[r]/=t}}return s}
/**
     * Returns the dot product. If {@link B} is an Array or Float64Array then an Array gets returned. If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */dot(t){if(t instanceof Matrix){let e=this;const[r,s]=e.shape,[i,n]=t.shape;if(s!==i)throw new Error(`A.dot(B): A is a ${e.shape.join(" ⨯ ")}-Matrix, B is a ${t.shape.join(" ⨯ ")}-Matrix:\n                A has ${s} cols and B ${i} rows.\n                Must be equal!`);return new Matrix(r,n,((r,i)=>{const o=e.row(r),a=t.values;let h=0;for(let t=0,e=i;t<s;++t,e+=n)h+=o[t]*a[e];return h}))}if(Matrix.isArray(t)){let e=this._rows;if(t.length!==e)throw new Error(`A.dot(B): A has ${e} cols and B has ${t.length} rows. Must be equal!`);let r=new Array(e);for(let s=0;s<e;++s)r[s]=neumair_sum(this.row(s).map((e=>e*t[s])));return r}throw new Error("B must be Matrix or Array")}
/**
     * Transposes the current matrix and returns the dot product with {@link B}.
     * If {@link B} is an Array or Float64Array then an Array gets returned.
     * If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */transDot(t){if(t instanceof Matrix){let e=this;const[r,s]=e.shape,[i,n]=t.shape;// transpose matrix
if(r!==i)throw new Error(`A.dot(B): A is a ${[s,r].join(" ⨯ ")}-Matrix, B is a ${t.shape.join(" ⨯ ")}-Matrix:\n                A has ${r} cols and B ${i} rows, which must be equal!`);
// let B = new Matrix(this._cols, this._rows, (row, col) => this.entry(col, row));
// this.values[row * this._cols + col];
return new Matrix(s,n,((i,o)=>{const a=e.values,h=t.values;let l=0;for(let t=0,e=i,_=o;t<r;++t,e+=s,_+=n)l+=a[e]*h[_];return l}))}if(Matrix.isArray(t)){let e=this._cols;if(t.length!==e)throw new Error(`A.dot(B): A has ${e} cols and B has ${t.length} rows. Must be equal!`);let r=new Array(e);for(let s=0;s<e;++s)r[s]=neumair_sum(this.col(s).map((e=>e*t[s])));return r}throw new Error("B must be Matrix or Array")}
/**
     * Returns the dot product with the transposed version of {@link B}.
     * If {@link B} is an Array or Float64Array then an Array gets returned.
     * If {@link B} is a Matrix then a Matrix gets returned.
     * @param {(Matrix|Array|Float64Array)} B the right side
     * @returns {(Matrix|Array)}
     */dotTrans(t){if(t instanceof Matrix){let e=this;const[r,s]=e.shape,[i,n]=t.shape;if(s!==n)throw new Error(`A.dot(B): A is a ${e.shape.join(" ⨯ ")}-Matrix, B is a ${[n,i].join(" ⨯ ")}-Matrix:\n                A has ${s} cols and B ${n} rows, which must be equal!`);return new Matrix(r,i,((r,i)=>{const n=e.row(r),o=t.row(i);let a=0;for(let t=0;t<s;++t)a+=n[t]*o[t];return a}))}if(Matrix.isArray(t)){let e=this._rows;if(t.length!==e)throw new Error(`A.dot(B): A has ${e} cols and B has ${t.length} rows. Must be equal!`);let r=new Array(e);for(let s=0;s<e;++s)r[s]=neumair_sum(this.row(s).map((e=>e*t[s])));return r}throw new Error("B must be Matrix or Array")}
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
     */concat(t,e="horizontal"){const r=this,[s,i]=r.shape,[n,o]=t.shape;if("horizontal"==e){if(s!=n)throw new Error(`A.concat(B, "horizontal"): A and B need same number of rows, A has ${s} rows, B has ${n} rows.`);const e=new Matrix(s,i+o,"zeros");return e.set_block(0,0,r),e.set_block(0,i,t),e}if("vertical"==e){if(i!=o)throw new Error(`A.concat(B, "vertical"): A and B need same number of columns, A has ${i} columns, B has ${o} columns.`);const e=new Matrix(s+n,i,"zeros");return e.set_block(0,0,r),e.set_block(s,0,t),e}if("diag"==e){const e=new Matrix(s+n,i+o,"zeros");return e.set_block(0,0,r),e.set_block(s,i,t),e}throw new Error(`type must be "horizontal" or "vertical", but type is ${e}!`)}
/**
     * Writes the entries of B in A at an offset position given by {@link offset_row} and {@link offset_col}.
     * @param {int} offset_row
     * @param {int} offset_col
     * @param {Matrix} B
     * @returns {Matrix}
     */set_block(t,e,r){const s=Math.min(this._rows-t,r.shape[0]),i=Math.min(this._cols-e,r.shape[1]);for(let n=0;n<s;++n)for(let s=0;s<i;++s)this.set_entry(n+t,s+e,r.entry(n,s));return this}
/**
     * Extracts the entries from the {@link start_row}<sup>th</sup> row to the {@link end_row}<sup>th</sup> row, the {@link start_col}<sup>th</sup> column to the {@link end_col}<sup>th</sup> column of the matrix.
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
     * A.get_block(1, 1); // [[5, 6], [8, 9]]
     * A.get_block(0, 0, 1, 1); // [[1]]
     * A.get_block(1, 1, 2, 2); // [[5]]
     * A.get_block(0, 0, 2, 2); // [[1, 2], [4, 5]]
     */get_block(t,e,r=null,s=null){const[i,n]=this.shape;if(s=s??n,(r=r??i)<=t||s<=e)throw new Error(`\n                end_row must be greater than start_row, and\n                end_col must be greater than start_col, but\n                end_row = ${r}, start_row = ${t}, end_col = ${s}, and start_col = ${e}!`);const o=new Matrix(r-t,s-e,"zeros");for(let i=t,n=0;i<r;++i,++n)for(let t=e,r=0;t<s;++t,++r)o.set_entry(n,r,this.entry(i,t));return o;
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
     * @private
     * @param {Function} f function takes 2 parameters, the value of the actual entry and a value given by the function {@link v}. The result of {@link f} gets writen to the Matrix.
     * @param {Function} v function takes 2 parameters for row and col, and returns a value witch should be applied to the colth entry of the rowth row of the matrix.
     */_apply_array(t,e){const r=this.values,[s,i]=this.shape;for(let n=0,o=0;o<s;++o)for(let s=0;s<i;++s,++n)r[n]=t(r[n],e(o,s));return this}_apply_rowwise_array(t,e){return this._apply_array(e,((e,r)=>t[r]))}_apply_colwise_array(t,e){const r=this.values,[s,i]=this.shape;for(let n=0,o=0;o<s;++o){const s=t[o];for(let t=0;t<i;++t,++n)r[n]=e(r[n],s)}return this}_apply(t,e){const r=this.values,[s,i]=this.shape;if(t instanceof Matrix){const n=t.values,[o,a]=t.shape;if(1===o){if(i!==a)throw new Error("cols !== value_cols");for(let t=0,o=0;o<s;++o)for(let s=0;s<i;++s,++t)r[t]=e(r[t],n[s])}else if(1===a){if(s!==o)throw new Error("rows !== value_rows");for(let t=0,o=0;o<s;++o){const s=n[o];for(let n=0;n<i;++n,++t)r[t]=e(r[t],s)}}else{if(s!=o||i!=a)throw new Error("error");for(let t=0,o=s*i;t<o;++t)r[t]=e(r[t],n[t])}}else if(Matrix.isArray(t))if(t.length===s)for(let n=0,o=0;o<s;++o){const s=t[o];for(let t=0;t<i;++t,++n)r[n]=e(r[n],s)}else{if(t.length!==i)throw new Error("error");for(let n=0,o=0;o<s;++o)for(let s=0;s<i;++s,++n)r[n]=e(r[n],t[s])}else// scalar value
for(let n=0,o=s*i;n<o;++n)r[n]=e(r[n],t);return this}
/**
     * Clones the Matrix.
     * @returns {Matrix}
     */clone(){let t=new Matrix;return t._rows=this._rows,t._cols=this._cols,t._data=this.values.slice(0),t}
/**
     * Entrywise multiplication with {@link value}.
     * @param {Matrix|Array|Number} value
     * @param {Object} [options]
     * @param {Boolean} [options.inline = false]  - If true, applies multiplication to the element, otherwise it creates first a copy and applies the multiplication on the copy.
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.mult(2); // [[2, 4], [6, 8]];
     * A.mult(B); // [[1, 4], [9, 16]];
     */mult(t,{inline:e=!1}={}){return(e?this:this.clone())._apply(t,((t,e)=>t*e))}
/**
     * Entrywise division with {@link value}.
     * @param {Matrix|Array|Number} value
     * @param {Object} [options]
     * @param {Boolean} [options.inline = false] - If true, applies division to the element, otherwise it creates first a copy and applies the division on the copy.
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.divide(2); // [[0.5, 1], [1.5, 2]];
     * A.divide(B); // [[1, 1], [1, 1]];
     */divide(t,{inline:e=!1}={}){return(e?this:this.clone())._apply(t,((t,e)=>t/e))}
/**
     * Entrywise addition with {@link value}.
     * @param {Matrix|Array|Number} value
     * @param {Object} [options]
     * @param {Boolean} [options.inline = false]  - If true, applies addition to the element, otherwise it creates first a copy and applies the addition on the copy.
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.add(2); // [[3, 4], [5, 6]];
     * A.add(B); // [[2, 4], [6, 8]];
     */add(t,{inline:e=!1}={}){return(e?this:this.clone())._apply(t,((t,e)=>t+e))}
/**
     * Entrywise subtraction with {@link value}.
     * @param {Matrix|Array|Number} value
     * @param {Object} [options]
     * @param {Boolean} [options.inline = false] - If true, applies subtraction to the element, otherwise it creates first a copy and applies the subtraction on the copy.
     * @returns {Matrix}
     * @example
     *
     * let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
     * let B = A.clone(); // B == A;
     *
     * A.sub(2); // [[-1, 0], [1, 2]];
     * A.sub(B); // [[0, 0], [0, 0]];
     */sub(t,{inline:e=!1}={}){return(e?this:this.clone())._apply(t,((t,e)=>t-e))}
/**
     * Returns the number of rows and columns of the Matrix.
     * @returns {Array} An Array in the form [rows, columns].
     */get shape(){return[this._rows,this._cols]}
/**
     * Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.
     * @param {Array} parameter - takes an Array in the form [rows, cols, value], where rows and cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row and col) which has to return a value for the colth entry of the rowth row.
     * @returns {Matrix}
     */set shape([t,e,r=(()=>0)]){this._rows=t,this._cols=e,this._data=new Float64Array(t*e);for(let s=0,i=0;i<t;++i)for(let t=0;t<e;++t,++s)this._data[s]=r(i,t);return this}
/**
     * Returns the Matrix as a Array of Float64Arrays.
     * @returns {Array<Float64Array>}
     */get to2dArray(){const t=[];for(const e of this.iterate_rows())t.push(e);return t}
/**
     * Returns the Matrix as a Array of Arrays.
     * @returns {Array<Array>}
     */get asArray(){const t=[];for(const e of this.iterate_rows())t.push(Array.from(e));return t}
/**
     * Returns the diagonal of the Matrix.
     * @returns {Float64Array}
     */get diag(){const t=this._rows,e=this._cols,r=Math.min(t,e);let s=new Float64Array(r);for(let t=0;t<r;++t)s[t]=this.entry(t,t);return s}
/**
     * Returns the mean of all entries of the Matrix.
     * @returns {Number}
     */get mean(){return this.sum/(this._rows*this._cols)}
/**
     * Returns the sum oof all entries of the Matrix.
     * @returns {Number}
     */get sum(){return neumair_sum(this.values)}
/**
     * Returns the entries of the Matrix.
     * @returns {Float64Array}
     */get values(){return this._data}
/**
     * Returns the mean of each row of the matrix.
     * @returns {Float64Array}
     */get meanRows(){const t=this.values,e=this._rows,r=this._cols,s=Float64Array.from({length:e});for(let i=0,n=0;n<e;++n){let e=0;for(let s=0;s<r;++s,++i)e+=t[i];s[n]=e/r}return s}
/** Returns the mean of each column of the matrix.
     * @returns {Float64Array}
     */get meanCols(){const t=this.values,e=this._rows,r=this._cols,s=Float64Array.from({length:r});for(let i=0;i<r;++i){let n=0;for(let s=i,o=0;o<e;++o,s+=r)n+=t[s];s[i]=n/e}return s}
/**
     * Solves the equation {@link A}x = {@link b} using the conjugate gradient method. Returns the result x.
     * @param {Matrix} A - Matrix
     * @param {Matrix} b - Matrix
     * @param {Randomizer} [randomizer=null]
     * @param {Number} [tol=1e-3]
     * @returns {Matrix}
     */static solve_CG(t,e,r,s=.001){null===r&&(r=new Randomizer);const i=t.shape[0],n=e.shape[1];let o=new Matrix(i,0);for(let a=0;a<n;++a){const n=Matrix.from(e.col(a)).T;let h=new Matrix(i,1,(()=>r.random)),l=n.sub(t.dot(h)),_=l.clone();do{const e=t.dot(_),r=l.transDot(l).entry(0,0)/_.transDot(e).entry(0,0);h=h.add(_.mult(r));const s=l.sub(e.mult(r)),i=s.transDot(s).entry(0,0)/l.transDot(l).entry(0,0);_=s.add(_.mult(i)),l=s}while(Math.abs(l.mean)>s);o=o.concat(h,"horizontal")}return o}
/**
     * Solves the equation {@link A}x = {@link b}. Returns the result x.
     * @param {Matrix} A - Matrix or LU Decomposition
     * @param {Matrix} b - Matrix
     * @returns {Matrix}
     */static solve(t,e){let{L:r,U:s}="L"in t&&"U"in t?t:Matrix.LU(t),i=r.shape[0],n=e.clone();
// forward
for(let t=0;t<i;++t){for(let e=0;e<t-1;++e)n.sub_entry(0,t,r.entry(t,e)*n.entry(1,e));n.set_entry(0,t,n.entry(0,t)/r.entry(t,t))}
// backward
for(let t=i-1;t>=0;--t){for(let e=i-1;e>t;--e)n.sub_entry(0,t,s.entry(t,e)*n.entry(0,e));n.set_entry(0,t,n.entry(0,t)/s.entry(t,t))}return n}
/**
     * {@link L}{@link U} decomposition of the Matrix {@link A}. Creates two matrices, so that the dot product LU equals A.
     * @param {Matrix} A
     * @returns {{L: Matrix, U: Matrix}} result - Returns the left triangle matrix {@link L} and the upper triangle matrix {@link U}.
     */static LU(t){const e=t.shape[0],r=new Matrix(e,e,"zeros"),s=new Matrix(e,e,"identity");for(let i=0;i<e;++i){for(let n=i;n<e;++n){let e=0;for(let t=0;t<i;++t)e+=r.entry(n,t)*s.entry(t,i);r.set_entry(n,i,t.entry(n,i)-e)}for(let n=i;n<e;++n){if(0===r.entry(i,i))return;let e=0;for(let t=0;t<i;++t)e+=r.entry(i,t)*s.entry(t,n);s.set_entry(i,n,(t.entry(i,n)-e)/r.entry(i,i))}}return{L:r,U:s}}
/**
     * Computes the determinante of {@link A}, by using the LU decomposition of {@link A}.
     * @param {Matrix} A
     * @returns {Number} det - Returns the determinate of the Matrix {@link A}.
     */static det(t){const e=t.shape[0],{L:r,U:s}=Matrix.LU(t),i=r.diag,n=s.diag;let o=i[0]*n[0];for(let t=1;t<e;++t)o*=i[t]*n[t];return o}
/**
     * Computes the {@link k} components of the SVD decomposition of the matrix {@link M}
     * @param {Matrix} M
     * @param {int} [k=2]
     * @returns {{U: Matrix, Sigma: Matrix, V: Matrix}}
     */static SVD(t,e=2){let r=t.transDot(t),s=t.dotTrans(t),{eigenvectors:i,eigenvalues:n}=simultaneous_poweriteration(r,e),{eigenvectors:o}=simultaneous_poweriteration(s,e);return{U:o,Sigma:n.map((t=>Math.sqrt(t))),V:i};
//Algorithm 1a: Householder reduction to bidiagonal form:
/* const [m, n] = A.shape;
        let U = new Matrix(m, n, (i, j) => i == j ? 1 : 0);
        console.log(U.to2dArray)
        let V = new Matrix(n, m, (i, j) => i == j ? 1 : 0);
        console.log(V.to2dArray)
        let B = Matrix.bidiagonal(A.clone(), U, V);
        console.log(U,V,B)
        return { U: U, "Sigma": B, V: V }; */}static isArray(t){return Array.isArray(t)||t instanceof Float64Array||t instanceof Float32Array}}
/**
 * @class
 * @memberof module:utils
 * @alias Randomizer
 */class Randomizer{
/**
     * Mersenne Twister random number generator.
     * @constructor
     * @param {Number} [_seed=new Date().getTime()] - The seed for the random number generator. If <code>_seed == null</code> then the actual time gets used as seed.
     * @see https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js
     */
constructor(t){return this._N=624,this._M=397,this._MATRIX_A=2567483615,this._UPPER_MASK=2147483648,this._LOWER_MASK=2147483647,this._mt=new Array(this._N),this._mti=this.N+1,this.seed=t||(new Date).getTime(),this}set seed(t){this._seed=t;let e=this._mt;for(e[0]=t>>>0,this._mti=1;this._mti<this._N;this._mti+=1){let t=this._mti,r=e[t-1]^e[t-1]>>>30;e[t]=(1812433253*((4294901760&r)>>>16)<<16)+1812433253*(65535&r)+t,e[t]>>>=0}}
/**
     * Returns the seed of the random number generator.
     * @returns {Number} - The seed.
     */get seed(){return this._seed}
/**
     * Returns a float between 0 and 1.
     * @returns {Number} - A random number between [0, 1]
     */get random(){return this.random_int*(1/4294967296)}
/**
     * Returns an integer between 0 and MAX_INTEGER.
     * @returns {Integer} - A random integer.
     */get random_int(){let t,e=new Array(0,this._MATRIX_A);if(this._mti>=this._N){let r,s=this._N-this._M,i=this._M-this._N;
/* if (this._mti == this._N + 1) {
                this.seed = 5489;
            } */for(r=0;r<s;++r)t=this._mt[r]&this._UPPER_MASK|this._mt[r+1]&this._LOWER_MASK,this._mt[r]=this._mt[r+this._M]^t>>>1^e[1&t];for(;r<this._N-1;++r)t=this._mt[r]&this._UPPER_MASK|this._mt[r+1]&this._LOWER_MASK,this._mt[r]=this._mt[r+i]^t>>>1^e[1&t];t=this._mt[this._N-1]&this._UPPER_MASK|this._mt[0]&this._LOWER_MASK,this._mt[this._N-1]=this._mt[this._M-1]^t>>>1^e[1&t],this._mti=0}return t=this._mt[this._mti+=1],t^=t>>>11,t^=t<<7&2636928640,t^=t<<15&4022730752,t^=t>>>18,t>>>0}gauss_random(){let t,e,r;if(null!=this._val)return t=this._val,this._val=null,t;do{t=2*this.random-1,e=2*this.random-1,r=t*t+e*e}while(!r||r>1);const s=Math.sqrt(-2*Math.log(r)/r);// cache this for next function call for efficiency
return this._val=e*s,t*s}
/**
     * Returns samples from an input Matrix or Array.
     * @param {Matrix|Array|Float64Array} A - The input Matrix or Array.
     * @param {Number} n - The number of samples.
     * @returns {Array} - A random selection form {@link A} of {@link n} samples.
     */choice(t,e){if(t instanceof Matrix){let r=t.shape[0];if(e>r)throw new Error("n bigger than A!");let s=new Array(e),i=linspace(0,r-1);for(let t=0,r=i.length;t<e;++t,--r){let e=this.random_int%r;s[t]=i.splice(e,1)[0]}return s.map((e=>t.row(e)))}if(Array.isArray(t)||t instanceof Float64Array){let r=t.length;if(e>r)throw new Error("n bigger than A!");let s=new Array(e),i=linspace(0,r-1);for(let t=0,r=i.length;t<e;++t,--r){let e=this.random_int%r;s[t]=i.splice(e,1)[0]}return s.map((e=>t[e]))}}
/**
     * @static
     * Returns samples from an input Matrix or Array.
     * @param {Matrix|Array|Float64Array} A - The input Matrix or Array.
     * @param {Number} n - The number of samples.
     * @param {Number} seed - The seed for the random number generator.
     * @returns {Array} - A random selection form {@link A} of {@link n} samples.
     */static choice(t,e,r=1212){return new Randomizer(r).choice(t,e);
/* let rows = A.shape[0];
        if (n > rows) {
            throw new Error("n bigger than A!");
        }
        let rand = new Randomizer(seed);
        let sample = new Array(n);
        let index_list = linspace(0, rows - 1);
        for (let i = 0, l = index_list.length; i < n; ++i, --l) {
            let random_index = rand.random_int % l;
            sample[i] = index_list.splice(random_index, 1)[0];
        }
        //return result;
        //return new Matrix(n, cols, (row, col) => A.entry(sample[row], col))
        return sample.map((d) => A.row(d)); */}}
/**
 * Returns maximum in Array {@link values}.
 * @memberof module:utils
 * @alias max
 * @param {Array} values 
 * @returns {Number}
 */function max(t){let e;for(const r of t)null!=r&&(e<r||void 0===e&&r>=r)&&(e=r);return e}
/**
 * Returns maximum in Array {@link values}.
 * @memberof module:utils
 * @alias min
 * @param {Array} values
 * @returns {Number}
 */function min(t){let e;for(const r of t)null!=r&&(e>r||void 0===e&&r<=r)&&(e=r);return e}
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
 * @class
 * @alias DR
 * @borrows DR#parameter as DR#para
 * @borrows DR#parameter as DR#p
 */class DR{
/**
     * Takes the default parameters and seals them, remembers the type of input {@link X}, and initializes the random number generator.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias DR
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed value for the random number generator.
     * @returns {DR}
     */
constructor(t,e,r){if(this._parameters=Object.assign(Object.seal(e),r),Array.isArray(t))this._type="array",this.X=Matrix.from(t);else{if(!(t instanceof Matrix))throw new Error("No valid type for X!");this._type="matrix",this.X=t}return[this._N,this._D]=this.X.shape,this._randomizer=new Randomizer(this._parameters.seed),this._is_initialized=!1,this}
/**
     * Set and get parameters
     * @param {String} [name = null] - Name of the parameter. If not given then returns all parameters as an Object.
     * @param {any} [value = null] - Value of the parameter to set. If <code>name</code> is set and <code>value</code> is not given, returns the value of the respective parameter.
     * @returns {DR|any|Object} 
     * On setting a parameter, this function returns the DR object. 
     * If <code>name</code> is set and <code>value == null</code> then return actual parameter value.
     * If <code>name</code> is not given, then returns all parameters as an Object.
     * 
     * @example
     * '''
     * const DR = new druid.TSNE(X, {d: 3}); // creates a new DR object, with parameter for <code>d</code> = 3.
     * DR.parameter("d"); // returns 3,
     * DR.parameter("d", 2); // sets parameter <code>d</code> to 2 and returns <code>DR</code>.
     * '''
     */parameter(t=null,e=null){if(null===t)return Object.assign({},this._parameters);if(!this._parameters.hasOwnProperty(t))throw new Error(`${t} is not a valid parameter!`);return null!==e?(this._parameters[t]=e,this._is_initialized=!1,this):this._parameters[t]}para(t=null,e=null){return this.parameter(t,e)}p(t=null,e=null){return this.parameter(t,e)}
/**
     * Computes the projection.
     * @returns {Matrix} the projection.
     */transform(){return this.check_init(),this.projection}
/**
     * Computes the projection.
     * @yields {Matrix|Number[][]} the intermediate steps of the projection.
     */*generator(){return this.transform()}
/**
     * If the respective DR method has an <code>init</code> function, call it before <code>transform</code>.
     * @returns {DR}
     */check_init(){return this._is_initialized||"function"!=typeof this.init||(this.init(),this._is_initialized=!0),this}
/**
     * @returns {Matrix|Number[][]} the projection in the type of input <code>X</code>.
     */get projection(){if(this.hasOwnProperty("Y"))return this.check_init(),"matrix"===this._type?this.Y:this.Y.to2dArray;throw new Error("The dataset is not transformed yet!")}
/**
     * Computes the projection.
     * @param  {...unknown} args - Arguments the transform method of the respective DR method takes.
     * @returns {Promise<Matrix|Number[][]>} the dimensionality reduced dataset.
     */async transform_async(...t){return this.transform(...t)}
/**
     * Computes the projection.
     * @static
     * @param  {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Matrix|Array} the dimensionality reduced dataset.
     */static transform(...t){return new this(...t).transform()}
/**
     * Computes the projection.
     * @static
     * @param  {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Promise} a promise yielding the dimensionality reduced dataset.
     */static async transform_async(...t){return this.transform(...t)}
/**
     * Computes the projection.
     * @static
     * @param  {...unknown} args - Takes the same arguments of the constructor of the respective DR method.
     * @returns {Generator} a generator yielding the intermediate steps of the dimensionality reduction method.
     */static*generator(...t){const e=new this(...t).generator();for(const t of e)yield t}}
/**
 * @class
 * @alias PCA
 * @augments DR
 */class PCA extends DR{
/**
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias PCA
     * @param {Matrix|Array<Array<Number>>} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @returns {PCA}
     */
constructor(t,e){return super(t,{d:2,seed:1212,eig_args:{}},e),this._parameters.eig_args.hasOwnProperty("seed")||(this._parameters.eig_args.seed=this._randomizer),this}
/**
     * Transforms the inputdata {@link X} to dimensionality {@link d}. If parameter {@link A} is given, then project {@link A} with the principal components of {@link X}.
     * @param {null|Matrix|Array} [A = null] - If given, the data to project.
     * @returns {Matrix|Array} - The projected data.
     */transform(t=null){const e=this.principal_components();if(null==t){const t=this.X;return this.Y=t.dot(e),this.projection}if(Array.isArray(t))return Matrix.from(t).dot(e).asArray;if(t instanceof Matrix)return t.dot(e);throw new Error("No valid type for A!")}
/**
     * Computes the {@link d} principal components of Matrix {@link X}.
     * @returns {Matrix}
     */principal_components(){if(this.V)return this.V;const{d:t,eig_args:e}=this._parameters,r=this.X,s=r.sub(r.meanCols),i=s.transDot(s),{eigenvectors:n}=simultaneous_poweriteration(i,t,e);return this.V=Matrix.from(n).transpose(),this.V}static principal_components(t,e){return new this(t,e).principal_components()}}
/**
 * @class
 * @alias MDS
 * @extends DR
 */class MDS extends DR{
/**
     * Classical MDS.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias MDS
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     */
constructor(t,e){return super(t,{d:2,metric:euclidean,seed:1212,eig_args:{}},e),this._parameters.eig_args.hasOwnProperty("seed")||(this._parameters.eig_args.seed=this._randomizer),this}
/**
     * Transforms the inputdata {@link X} to dimensionality {@link d}.
     * @returns {Matrix|Array}
     */transform(){const t=this.X,e=t.shape[0],{d:r,metric:s,eig_args:i}=this._parameters,n="precomputed"===s?t:distance_matrix(t,s),o=n.meanCols,a=n.meanRows,h=n.mean;this._d_X=n;const l=new Matrix(e,e,((t,e)=>n.entry(t,e)-o[t]-a[e]+h)),{eigenvectors:_}=simultaneous_poweriteration(l,r,i);return this.Y=Matrix.from(_).transpose(),this.projection}
/**
     * @returns {Number} - the stress of the projection.
     */stress(){const t=this.X.shape[0],e=this.Y,r=this._d_X,s=new Matrix;s.shape=[t,t,(t,r)=>t<r?euclidean(e.row(t),e.row(r)):s.entry(r,t)];let i=0,n=0;for(let e=0;e<t;++e)for(let o=e+1;o<t;++o)i+=Math.pow(r.entry(e,o)-s.entry(e,o),2),n+=Math.pow(r.entry(e,o),2);return Math.sqrt(i/n)}}
/**
 * @class
 * @alias ISOMAP
 * @extends DR
 */class ISOMAP extends DR{
/**
     * Isometric feature mapping (ISOMAP).
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias ISOMAP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} parameters.neighbors - the number of neighbors {@link ISOMAP} should use to project the data.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://doi.org/10.1126/science.290.5500.2319}
     */
constructor(t,e){return super(t,{neighbors:void 0,d:2,metric:euclidean,seed:1212,eig_args:{}},e),this.parameter("neighbors",Math.min(this._parameters.neighbors??Math.max(Math.floor(this.X.shape[0]/10),2),this._N-1)),this._parameters.eig_args.hasOwnProperty("seed")||(this._parameters.eig_args.seed=this._randomizer),this}
/**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */transform(){this.check_init();const t=this.X,e=this._N,{d:r,metric:s,eig_args:i,neighbors:n}=this._parameters,o=new Matrix;o.shape=[e,e,(e,r)=>e<=r?s(t.row(e),t.row(r)):o.entry(r,e)];const a=[];for(let t=0;t<e;++t){const r=[];for(let s=0;s<e;++s)r.push({index:s,distance:o.entry(t,s)});const s=new Heap(r,(t=>t.distance),"min");a.push(s.toArray().slice(1,n+1))}
/*D = dijkstra(kNearestNeighbors);*/
// compute shortest paths
// TODO: make extern
/** @see {@link https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm} */const h=new Matrix(e,e,((t,e)=>{const r=a[t].find((t=>t.index===e));return r?r.distance:1/0}));for(let t=0;t<e;++t)for(let r=0;r<e;++r){let s=h.entry(t,r);for(let i=0;i<e;++i)s=Math.min(s,h.entry(t,i)+h.entry(i,r));h.set_entry(t,r,s)}let l=new Float64Array(e),_=new Float64Array(e),c=0;const u=new Matrix(e,e,((t,e)=>{let r=h.entry(t,e);return r=r===1/0?0:r,l[t]+=r,_[e]+=r,c+=r,r}));l=l.map((t=>t/e)),_=_.map((t=>t/e)),c/=e**2;const d=new Matrix(e,e,((t,e)=>u.entry(t,e)-l[t]-_[e]+c)),{eigenvectors:m}=simultaneous_poweriteration(d,r,i);
// compute d eigenvectors
// return embedding
return this.Y=Matrix.from(m).transpose(),this.projection}}
/**
 * @class
 * @alias FASTMAP
 * @extends DR
 */class FASTMAP extends DR{
/**
     * FastMap: a fast algorithm for indexing, data-mining and visualization of traditional and multimedia datasets
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias FASTMAP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the dimensionality of the projection.
     * @returns {FASTMAP}
     * @see {@link https://doi.org/10.1145/223784.223812}
     */
constructor(t,e){return super(t,{d:2,metric:euclidean,seed:1212},e),this}
/**
     * Chooses two points which are the most distant in the actual projection.
     * @private
     * @param {Function} dist
     * @returns {Array} An array consisting of first index, second index, and distance between the two points.
     */_choose_distant_objects(t){const e=this.X.shape[0];let r=this._randomizer.random_int%e-1,s=null,i=-1/0;for(let n=0;n<e;++n){const e=t(r,n);e>i&&(i=e,s=n)}i=-1/0;for(let n=0;n<e;++n){const e=t(s,n);e>i&&(i=e,r=n)}return[r,s,i]}
/**
     * Computes the projection.
     * @returns {Matrix} The {@link d}-dimensional projection of the data matrix {@link X}.
     */transform(){const t=this.X,e=t.shape[0],{d:r,metric:s}=this._parameters,i=new Matrix(e,r,0);let dist=(e,r)=>s(t.row(e),t.row(r));for(let t=0;t<r;++t){let r=dist;
// choose pivot objects
const[s,n,o]=this._choose_distant_objects(dist);if(0!==o){
// project the objects on the line (O_a, O_b)
for(let r=0;r<e;++r){const e=(dist(s,r)**2+o**2-dist(n,r)**2)/(2*o);i.set_entry(r,t,e)}
// consider the projections of the objects on a
// hyperplane perpendicluar to the line (a, b);
// the distance function D'() between two
// projections is given by Eq.4
dist=(e,s)=>Math.sqrt(r(e,s)**2-(i.entry(e,t)-i.entry(s,t))**2)}}
// return embedding.
return this.Y=i,this.projection}}
/**
 * @class
 * @alias LDA
 * @extends DR
 */class LDA extends DR{
/**
     * Linear Discriminant Analysis.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LDA
     * @param {Matrix} X - The high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Array} parameters.labels - The labels / classes for each data point.
     * @param {number} [parameters.d = 2] - The dimensionality of the projection.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x}
     */
constructor(t,e){return super(t,{labels:null,d:2,seed:1212,eig_args:{}},e),this._parameters.eig_args.hasOwnProperty("seed")||(this._parameters.eig_args.seed=this._randomizer),this}
/**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */transform(){const t=this.X,[e,r]=t.shape,{d:s,labels:i,eig_args:n}=this._parameters;if(null===i||i.length!=e)throw new Error("LDA needs parameter label to every datapoint to work!");const o={};let a=0;i.forEach(((e,r)=>{e in o?(o[e].count++,o[e].rows.push(t.row(r))):o[e]={id:a++,count:1,rows:[t.row(r)]}}));
// create X_mean and vector means;
const h=t.mean,l=new Matrix(a,r);for(const t in o){const e=Matrix.from(o[t].rows).meanCols;for(let s=0;s<r;++s)l.set_entry(o[t].id,s,e[s])}
// scatter_between
let _=new Matrix(r,r);for(const t in o){const e=l.row(o[t].id),s=new Matrix(r,1,(t=>e[t]-h)),i=o[t].count;_=_.add(s.dotTrans(s).mult(i))}
// scatter_within
let c=new Matrix(r,r);for(const t in o){const e=l.row(o[t].id),s=new Matrix(r,1,(t=>e[t])),i=o[t].rows;for(let e=0,n=o[t].count;e<n;++e){const t=new Matrix(r,1,((t,r)=>i[e][t]-s.entry(t,0)));c=c.add(t.dotTrans(t))}}let{eigenvectors:u}=simultaneous_poweriteration(c.inverse().dot(_),s,n);
// return embedding
return u=Matrix.from(u).transpose(),this.Y=t.dot(u),this.projection}}
/**
 * @class
 * @alias LLE
 * @extends DR
 */class LLE extends DR{
/**
     * Locally Linear Embedding.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LLE
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} neighbors - the label / class of each data point.
     * @param {Number} [d = 2] - the dimensionality of the projection.
     * @param {Function} [metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [seed = 1212] - the dimensionality of the projection.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://doi.org/10.1126/science.290.5500.2323}
     */
constructor(t,e){return super(t,{neighbors:void 0,d:2,metric:euclidean,seed:1212,eig_args:{}},e),this.parameter("neighbors",Math.min(e.neighbors??Math.max(Math.floor(this._N/10),2),this._N-1)),this._parameters.eig_args.hasOwnProperty("seed")||(this._parameters.eig_args.seed=this._randomizer),this}
/**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */transform(){const t=this.X,e=this._N,r=this._D,{neighbors:s,d:i,eig_args:n,metric:o}=this._parameters,a=k_nearest_neighbors(t,s,o),h=new Matrix(s,1,1),l=new Matrix(e,e);for(let i=0;i<e;++i){const e=a[i],n=new Matrix(s,r,((r,s)=>t.entry(e[r].j,s)-t.entry(i,s))),o=n.dotTrans(n);if(s>r){const t=neumair_sum(o.diag)/1e3;for(let e=0;e<s;++e)o.add_entry(e,e,t)}
// reconstruct;
let _=Matrix.solve_CG(o,h,this._randomizer);_=_.divide(_.sum);for(let t=0;t<s;++t)l.set_entry(i,e[t].j,_.entry(t,0))}
// comp embedding
const _=new Matrix(e,e,"identity").sub(l),c=_.transDot(_),{eigenvectors:u}=simultaneous_poweriteration(c.T.inverse(),i+1,n);
// return embedding
return this.Y=Matrix.from(u.slice(1,1+i)).T,this.projection}}
/**
 * @class
 * @alias LTSA
 * @extends DR
 */class LTSA extends DR{
/**
     * Local Tangent Space Alignment
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LTSA
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} parameters.neighbors - the number of neighbors {@link LTSA} should use to project the data.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @param {Number} [parameters.eig_args] - Parameters for the eigendecomposition algorithm.
     * @see {@link https://epubs.siam.org/doi/abs/10.1137/S1064827502419154}
     */
constructor(t,e){if(super(t,{neighbors:void 0,d:2,metric:euclidean,seed:1212,eig_args:{}},e),this.parameter("neighbors",Math.min(e.neighbors??Math.max(Math.floor(this._N/10),2),this._N-1)),this._parameters.eig_args.hasOwnProperty("seed")||(this._parameters.eig_args.seed=this._randomizer),this._D<=this.parameter("d"))throw new Error(`Dimensionality of X (D = ${this._D}) must be greater than the required dimensionality of the result (d = ${this.parameter("d")})!`);return this}
/**
     * Transforms the inputdata {@link X} to dimenionality {@link d}.
     */transform(){const t=this.X,[e,r]=t.shape,{d:s,neighbors:i,metric:n,eig_args:o}=this._parameters,a=k_nearest_neighbors(t,i,n),h=new Matrix(r,r,"center"),l=new Matrix(e,e,0);for(let r=0;r<e;++r){
// 1.2 compute the d largest eigenvectors of the correlation matrix
const e=[r,...a[r].map((t=>t.j))];let n=Matrix.from(e.map((e=>t.row(e))));
// center X_i
n=n.dot(h);
// correlation matrix
const _=n.dotTrans(n),{eigenvectors:c}=simultaneous_poweriteration(_,s,o),u=Matrix.from(c),d=u.transDot(u).add(1/Math.sqrt(i+1));for(let t=0;t<i+1;++t)for(let r=0;r<i+1;++r)l.add_entry(e[t],e[r],d.entry(t,r)-(t===r?1:0))}
// 3. Aligning global coordinates
const{eigenvectors:_}=simultaneous_poweriteration(l,s+1,o);
// return embedding
return this.Y=Matrix.from(_.slice(1)).transpose(),this.projection}}
/**
 * @class
 * @alias TSNE
 * @extends DR
 */class TSNE extends DR{
/**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TSNE
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.perplexity = 50] - perplexity.
     * @param {Number} [parameters.epsilon = 10] - learning parameter.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [parameters.metric = euclidean_squared] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TSNE}
     */
constructor(t,e){return super(t,{perplexity:50,epsilon:10,d:2,metric:euclidean_squared,seed:1212},e),[this._N,this._D]=this.X.shape,this._iter=0,this.Y=new Matrix(this._N,this.parameter("d"),(()=>1e-4*this._randomizer.gauss_random())),this}
/**
     *
     * @returns {TSNE}
     */init(){
// init
const t=Math.log(this.parameter("perplexity")),e=this._N,r=this._D,{metric:s}=this._parameters,i=this.X;let n;if("precomputed"==s)n=druid.Matrix.from(i);else{n=new Matrix(e,e);for(let t=0;t<e;++t){const r=i.row(t);for(let o=t+1;o<e;++o){const e=s(r,i.row(o));n.set_entry(t,o,e),n.set_entry(o,t,e)}}}const o=new Matrix(e,e,0);this._ystep=new Matrix(e,r,0),this._gains=new Matrix(e,r,1);for(let r=0;r<e;++r){const s=n.row(r),i=o.row(r);let a,h=-1/0,l=1/0,_=1,c=50,u=!1;for(;!u&&c--;){
// compute entropy and kernel row with beta precision
a=0;let n=0;for(let t=0;t<e;++t){const e=s[t],o=r!==t?Math.exp(-e*_):0;n+=e*o,i[t]=o,a+=o}
// compute entropy
const o=a>0?Math.log(a)+_*n/a:0;o>t?(h=_,_=l===1/0?2*_:(_+l)/2):(l=_,_=h===-1/0?_/2:(_+h)/2),u=Math.abs(o-t)<1e-4}
// normalize p
for(let t=0;t<e;++t)i[t]/=a}
// compute probabilities
const a=2*e;for(let t=0;t<e;++t)for(let r=t;r<e;++r){const e=Math.max((o.entry(t,r)+o.entry(r,t))/a,1e-100);o.set_entry(t,r,e),o.set_entry(r,t,e)}return this._P=o,this}
/**
     *
     * @param {Number} [iterations=500] - Number of iterations.
     * @returns {Matrix|Number[][]} the projection.
     */transform(t=500){this.check_init();for(let e=0;e<t;++e)this.next();return this.projection}
/**
     *
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Number[][]} - the projection.
     */*generator(t=500){this.check_init();for(let e=0;e<t;++e)this.next(),yield this.projection;return this.projection}
/**
     * performs a optimization step
     * @private
     * @returns {Matrix}
     */next(){const t=++this._iter,e=this._P,r=this._ystep,s=this._gains,i=this._N,{d:n,epsilon:o}=this._parameters;let a=this.Y;
//calc cost gradient;
const h=t<100?4:1,l=new Matrix(i,i,"zeros");
// compute Q dist (unnormalized)
let _=0;for(let t=0;t<i;++t)for(let e=t+1;e<i;++e){let r=0;for(let s=0;s<n;++s){const i=a.entry(t,s)-a.entry(e,s);r+=i*i}const s=1/(1+r);l.set_entry(t,e,s),l.set_entry(e,t,s),_+=2*s}
// normalize Q dist
const c=new Matrix(i,i,0);for(let t=0;t<i;++t)for(let e=t+1;e<i;++e){const r=Math.max(l.entry(t,e)/_,1e-100);c.set_entry(t,e,r),c.set_entry(e,t,r)}const u=new Matrix(i,n,"zeros");for(let t=0;t<i;++t)for(let r=0;r<i;++r){const s=4*(h*e.entry(t,r)-c.entry(t,r))*l.entry(t,r);for(let e=0;e<n;++e)u.add_entry(t,e,s*(a.entry(t,e)-a.entry(r,e)))}
// perform gradient step
let d=new Float64Array(n);for(let e=0;e<i;++e)for(let i=0;i<n;++i){const n=u.entry(e,i),h=r.entry(e,i),l=s.entry(e,i);let _=Math.sign(n)===Math.sign(h)?.8*l:l+.2;_<.01&&(_=.01),s.set_entry(e,i,_);const c=(t<250?.5:.8)*h-o*_*n;r.set_entry(e,i,c),a.add_entry(e,i,c),d[i]+=a.entry(e,i)}for(let t=0;t<i;++t)for(let e=0;e<n;++e)a.sub_entry(t,e,d[e]/i);return this.Y}}
/**
 *
 * @memberof module:optimization
 * @alias powell
 * @param {Function} f
 * @param {Array} x0
 * @param {Number} [max_iter = 300]
 * @returns {Array}
 * @see http://optimization-js.github.io/optimization-js/optimization.js.html#line438
 */function powell(t,e,r=300){const s=e.length;let i=.001,n=1e4,o=e.slice(),a=t(o),h=!1;for(;r-- >=0&&!h;){h=!0;for(let e=0;e<s;++e){o[e]+=1e-6;let r=t(o);o[e]-=1e-6;let s=(r-a)/1e-6;Math.abs(s)>.01&&(h=!1),o[e]-=i*s,a=t(o)}i*=n>=a?1.05:.4,n=a}return o}
/**
 * @class
 * @alias UMAP
 * @extends DR
 */class UMAP extends DR{
/**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias UMAP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.n_neighbors = 15] - size of the local neighborhood.
     * @param {Number} [parameters.local_connectivity = 1] - number of nearest neighbors connected in the local neighborhood.
     * @param {Number} [parameters.min_dist = 1] - controls how tightly points get packed together.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points in the high-dimensional space.
     * @param {Number} [parameters._spread = 1] - The effective scale of embedded points. (In combination with {@link parameters.min_dist})
     * @param {Number} [parameters._set_op_mix_ratio = 1] - Interpolate between union and intersection.
     * @param {Number} [parameters._repulsion_strength = 1]  - Weighting applied to negative samples.
     * @param {Number} [parameters._negative_sample_rate = 5] - The number of negative samples per positive sample.
     * @param {Number} [parameters._n_epochs = 350] - The number of training epochs.
     * @param {Number} [parameter._initial_alpha = 1] - The initial learning rate for the optimization.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {UMAP}
     */
constructor(t,e){
/* let n_neighbors = Math.min(this._N - 1, parameters.n_neighbors);
        this.parameter("n_neighbors", n_neighbors);
        this.parameter("local_connectivity", Math.min(this.parameter("local_connectivity"), n_neighbors - 1)); */
if(super(t,{n_neighbors:15,local_connectivity:1,min_dist:1,d:2,metric:euclidean,seed:1212,_spread:1,_set_op_mix_ratio:1,_repulsion_strength:1,_negative_sample_rate:5,_n_epochs:350,_initial_alpha:1},e),[this._N,this._D]=this.X.shape,this.parameter("n_neighbors")>this._N)throw new Error(`Parameter n_neighbors (=${this.parameter("n_neighbors")}) needs to be smaller than dataset size (N=${this._N})!`);if(this.parameter("local_connectivity")>this.parameter("n_neighbors"))throw new Error(`Parameter local_connectivity (=${this.parameter("local_connectivity")}) needs to be smaller than parameter n_neighbors (=${this.parameter("n_neighbors")})`);this._iter=0;const r=this._randomizer;return this.Y=new Matrix(this._N,this.parameter("d"),(()=>r.random)),this}
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
     */_compute_membership_strengths(t,e,r){for(let s=0,i=t.length;s<i;++s){const i=r[s],n=t[s];for(let t=0,r=n.length;t<r;++t){const r=n[t].value-i;n[t].value=r>0?Math.exp(-r/e[s]):1}}return t}
/**
     * @private
     * @param {KNN|BallTree} knn
     * @param {Number} k
     * @returns {Object}
     */_smooth_knn_dist(t,e){const r=1e-5,s=.001,{local_connectivity:i,metric:n}=this._parameters,o=Math.log2(e),a=[],h=[],l=this.X,_=l.shape[0],c=[];if("precomputed"===n)for(let r=0;r<_;++r)c.push(t.search(r,e).reverse());else for(const r of l)c.push(t.search(r,e).raw_data().reverse());const u=Math.floor(i),d=i-u;for(let t=0;t<_;++t){let n=0,l=1/0,_=1,m=0;const f=c[t],p=f.filter((t=>t.value>0)),w=p.length;w>=i?u>0?(m=p[u-1].value,d>r&&(m+=d*(p[u].value-p[u-1].value))):m=d*p[0].value:w>0&&(m=p[w-1].value);for(let t=0;t<64;++t){let t=0;for(let r=0;r<e;++r){const e=f[r].value-m;t+=e>0?Math.exp(-e/_):1}if(Math.abs(t-o)<r)break;t>o?[l,_]=[_,(n+l)/2]:[n,_]=l===1/0?[_,2*_]:[_,(n+l)/2]}
//let mean_d = null;
if(m>0){const t=f.reduce(((t,e)=>t+e.value),0)/f.length;_<s*t&&(_=s*t)}else{const t=c.reduce(((t,e)=>t+e.reduce(((t,e)=>t+e.value),0)/e.length));_<s*t&&(_=s*t)}a[t]=m,h[t]=_}return{distances:c,sigmas:h,rhos:a}}
/**
     * @private
     * @param {Matrix} X
     * @param {Number} n_neighbors
     * @returns {Matrix}
     */_fuzzy_simplicial_set(t,e){const r=t.shape[0],{metric:s,_set_op_mix_ratio:i}=this._parameters,n="precomputed"===s?new KNN(t,"precomputed"):new BallTree(t.to2dArray,s);let{distances:o,sigmas:a,rhos:h}=this._smooth_knn_dist(n,e);o=this._compute_membership_strengths(o,a,h);const l=new Matrix(r,r,"zeros");for(let t=0;t<r;++t){const e=o[t];for(let r=0;r<e.length;++r)l.set_entry(t,e[r].element.index,e[r].value)}const _=l.T,c=l.mult(_);return l.add(_).sub(c).mult(i).add(c.mult(1-i))}
/**
     * @private
     * @param {Number} n_epochs
     * @returns {Array}
     */_make_epochs_per_sample(t){const e=this._weights,r=new Float32Array(e.length).fill(-1),s=t/max(e);return e.forEach(((e,i)=>{const n=e*s;n>0&&(r[i]=Math.round(t/n))})),r}
/**
     * @private
     * @param {Matrix} graph
     * @returns {Object}
     */_tocoo(t){const e=[],r=[],s=[],[i,n]=t.shape;for(let o=0;o<i;++o)for(let i=0;i<n;++i){const n=t.entry(o,i);0!==n&&(e.push(o),r.push(i),s.push(n))}return{rows:e,cols:r,data:s}}
/**
     * Computes all necessary
     * @returns {UMAP}
     */init(){const{_spread:t,min_dist:e,n_neighbors:r,_n_epochs:s,_negative_sample_rate:i}=this._parameters,[n,o]=this._find_ab_params(t,e);this._a=n,this._b=o,this._graph=this._fuzzy_simplicial_set(this.X,r);const{rows:a,cols:h,data:l}=this._tocoo(this._graph);return this._head=a,this._tail=h,this._weights=l,this._epochs_per_sample=this._make_epochs_per_sample(s),this._epochs_per_negative_sample=this._epochs_per_sample.map((t=>t*i)),this._epoch_of_next_sample=this._epochs_per_sample.slice(),this._epoch_of_next_negative_sample=this._epochs_per_negative_sample.slice(),this}graph(){return this.check_init(),{cols:this._head,rows:this._tail,weights:this._weights}}
/**
     *
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */transform(t=350){this.parameter("_n_epochs")!=t&&(this.parameter("_n_epochs",t),this.init()),this.check_init();for(let e=0;e<t;++e)this.next();return this.projection}
/**
     *
     * @param {Number} [iterations=350] - number of iterations.
     * @returns {Matrix|Array}
     */*generator(t=350){this.parameter("_n_epochs")!=t&&(this.parameter("_n_epochs",t),this.init()),this.check_init();for(let e=0;e<t;++e)this.next(),yield this.projection;return this.projection}
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
     */_optimize_layout(t,e,r,s){const i=this._randomizer,{_repulsion_strength:n,d:o}=this._parameters,{_alpha:a,_a:h,_b:l,_epochs_per_sample:_,_epochs_per_negative_sample:c,_epoch_of_next_negative_sample:u,_epoch_of_next_sample:d,_clip:m}=this,f=s.length;for(let p=0,w=_.length;p<w;++p)if(d[p]<=this._iter){const w=r[p],g=s[p],y=t.row(w),M=e.row(g),x=euclidean_squared(y,M);if(x>0){const t=-2*h*l*Math.pow(x,l-1)/(h*Math.pow(x,l)+1);for(let e=0;e<o;++e){const r=m(t*(y[e]-M[e]))*a;y[e]+=r,M[e]-=r}}d[p]+=_[p];const A=(this._iter-u[p])/c[p];for(let t=0;t<A;++t){const t=i.random_int%f,r=e.row(s[t]),_=euclidean_squared(y,r);if(_>0){const t=2*n*l/((.01+_)*(h*Math.pow(_,l)+1));for(let e=0;e<o;++e){const s=m(t*(y[e]-r[e]))*a;y[e]+=s,r[e]-=s}}else if(w===t)continue}u[p]+=A*c[p]}return t}
/**
     * @private
     * @returns {Matrix}
     */next(){const t=++this._iter,e=this.Y,{_initial_alpha:r,_n_epochs:s}=this._parameters;return this._alpha=r*(1-t/s),this.Y=this._optimize_layout(e,e,this._head,this._tail),this.Y}}
/**
 * @class
 * @alias TriMap
 * @extends DR
 */class TriMap extends DR{
/**
     *
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TriMap
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.weight_adj = 500] - scaling factor.
     * @param {Number} [parameters.c = 5] - number of triplets multiplier.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Number} [parameters.tol = 1e-8] -
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TriMap}
     * @see {@link https://arxiv.org/pdf/1910.00204v1.pdf}
     * @see {@link https://github.com/eamid/trimap}
     */
constructor(t,e){return super(t,{weight_adj:500,c:5,d:2,metric:euclidean,tol:1e-8,seed:1212},e),this}
/**
     *
     * @param {Matrix} [pca = null] - Initial Embedding (if null then PCA gets used).
     * @param {KNN} [knn = null] - KNN Object (if null then BallTree gets used).
     */init(t=null,e=null){const r=this.X,s=r.shape[0],{c:i,d:n,metric:o,seed:a}=this._parameters;this.n_inliers=2*i,this.n_outliers=1*i,this.n_random=1*i,this.Y=t||new PCA(r,{d:n,seed:a}).transform(),this.knn=e||new BallTree(r.to2dArray,o);const{triplets:h,weights:l}=this._generate_triplets(this.n_inliers,this.n_outliers,this.n_random);return this.triplets=h,this.weights=l,this.lr=1e3*s/h.shape[0],this.C=1/0,this.vel=new Matrix(s,n,0),this.gain=new Matrix(s,n,1),this}
/**
     * Generates {@link n_inliers} x {@link n_outliers} x {@link n_random} triplets.
     * @param {Number} n_inliers
     * @param {Number} n_outliers
     * @param {Number} n_random
     */_generate_triplets(t,e,r){const{metric:s,weight_adj:i}=this._parameters,n=this.X,o=n.shape[0],a=this.knn,h=Math.min(t+20,o),l=new Matrix(o,h),_=new Matrix(o,h);for(let t=0;t<o;++t)a.search(n.row(t),h+1).raw_data().filter((t=>0!=t.value)).sort(((t,e)=>t.value-e.value)).forEach(((e,r)=>{l.set_entry(t,r,e.element.index),_.set_entry(t,r,e.value)}));
// scale parameter
const c=new Float64Array(o);for(let t=0;t<o;++t)c[t]=Math.max((_.entry(t,3)+_.entry(t,4)+_.entry(t,5)+_.entry(t,6))/4,1e-10);const u=this._find_p(_,c,l);let d=this._sample_knn_triplets(u,l,t,e),m=d.shape[0];const f=new Float64Array(m);for(let t=0;t<m;++t){const e=d.entry(t,0),r=d.entry(t,2);f[t]=s(n.row(e),n.row(r))}let p=this._find_weights(d,u,l,f,c);if(r>0){const{random_triplets:t,random_weights:e}=this._sample_random_triplets(n,r,c);d=d.concat(t,"vertical"),p=Float64Array.from([...p,...e])}m=d.shape[0];let w=-1/0;for(let t=0;t<m;++t)isNaN(p[t])&&(p[t]=0),w<p[t]&&(w=p[t]);let g=-1/0;for(let t=0;t<m;++t)p[t]/=w,p[t]+=1e-4,p[t]=Math.log(1+i*p[t]),g<p[t]&&(g=p[t]);for(let t=0;t<m;++t)p[t]/=g;return{triplets:d,weights:p}}
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
     */_sample_knn_triplets(t,e,r,s){const i=e.shape[0],n=new Matrix(i*r*s,3);for(let o=0;o<i;++o){let a=o*r*s;const h=this.__argsort(t.row(o));for(let t=0;t<r;++t){let r=t*s;const l=e.entry(o,h[t]),_=this._rejection_sample(s,i,h.slice(0,t+1));for(let t=0;t<s;++t){const e=a+r+t,s=_[t];n.set_entry(e,0,o),n.set_entry(e,1,l),n.set_entry(e,2,s)}}}return n}
/**
     * Should do the same as np.argsort()
     * @private
     * @param {Array} A
     */__argsort(t){return linspace(0,t.length-1).sort(((e,r)=>t[r]-t[e]))}
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
     */_sample_random_triplets(t,e,r){const s=this.parameter("metric"),i=this._randomizer,n=t.shape[0],o=new Matrix(n*e,3),a=new Float64Array(n*e);for(let h=0;h<n;++h){const l=h*e,_=[...linspace(0,h-1),...linspace(h+1,n-1)];for(let n=0;n<e;++n){let[e,c]=i.choice(_,2),u=Math.exp(-(s(t.row(h),t.row(e))**2)/(r[h]*r[e]));u<1e-20&&(u=1e-20);let d=Math.exp(-(s(t.row(h),t.row(c))**2)/(r[h]*r[c]));d<1e-20&&(d=1e-20),u<d&&([e,c]=[c,e],[u,d]=[d,u]);const m=l+n;o.set_entry(m,0,h),o.set_entry(m,1,e),o.set_entry(m,2,c),a[m]=u/d}}return{random_triplets:o,random_weights:a}}
/**
     * Computes the gradient for updating the embedding.
     * @param {Matrix} Y - The embedding
     */_grad(t){const e=this.n_inliers,r=this.n_outliers,s=this.triplets,i=this.weights,[n,o]=t.shape,a=s.shape[0],h=new Matrix(n,o,0);let l=new Float64Array(o),_=new Float64Array(o),c=1,u=1,d=0,m=0;const f=n*e*r;for(let e=0;e<a;++e){const[n,a,p]=s.row(e);
// update y_ij, y_ik, d_ij, d_ik
if(e%r==0||e>=f){c=1,u=1;for(let e=0;e<o;++e){const r=t.entry(n,e),s=t.entry(a,e),i=t.entry(p,e);l[e]=r-s,_[e]=r-i,c+=l[e]**2,u+=_[e]**2}
// update y_ik and d_ik only
}else{u=1;for(let e=0;e<o;++e){const r=t.entry(n,e),s=t.entry(p,e);_[e]=r-s,u+=_[e]**2}}c>u&&++d,m+=i[e]/(1+u/c);const w=(i[e]/(c+u))**2;for(let t=0;t<o;++t){const e=l[t]*u*w,r=_[t]*c*w;h.add_entry(n,t,e-r),h.sub_entry(a,t,e),h.add_entry(p,t,r)}}return{grad:h,loss:m,n_viol:d}}
/**
     *
     * @param {Number} max_iteration
     */transform(t=400){this.check_init();for(let e=0;e<t;++e)this._next(e);return this.projection}
/**
     * @param {Number} max_iteration
     * @yields {Matrix}
     * @returns {Matrix}
     */*generator(t=800){this.check_init();for(let e=0;e<t;++e)this._next(e),yield this.projection;return this.projection}
/**
     * Does the iteration step.
     * @private
     * @param {Number} iter
     */_next(t){const e=t>150?.5:.3,r=this.C,s=this.vel,i=this.Y.add(s.mult(e)),{grad:n,loss:o,n_viol:a}=this._grad(i);return this.C=o,this.Y=this._update_embedding(i,t,n),this.lr*=r>o+this._parameters.tol?1.01:.9,this.Y}
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
constructor(t,e="complete",r=euclidean){if(this._id=0,this._matrix=t instanceof Matrix?t:Matrix.from(t),this._metric=r,this._linkage=e,"precomputed"===r&&this._matrix.shape[0]!==this._matrix.shape[1])throw new Error("If metric is 'precomputed', then matrix has to be square!");return this.init(),this.root=this.do(),this}
/**
     *
     * @param {Number} value - value where to cut the tree.
     * @param {("distance"|"depth")} [type = "distance"] - type of value.
     * @returns {Array<Array>} - Array of clusters with the indices of the rows in given {@link matrix}.
     */get_clusters(t,e="distance"){let r,s=[];switch(e){case"distance":r=t=>t.dist;break;case"depth":r=t=>t.depth;break;default:throw new Error("invalid type")}return this._traverse(this.root,r,t,s),s}
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
     */do(){const t=this._n,e=this._d_min,r=this._distance_matrix,s=this._clusters,i=this._c_size,n=this._linkage;let o=null;for(let a=0,h=t-1;a<h;++a){let a=0;for(let s=0;s<t;++s){let i=r.entry(s,e[s]);for(let n=s+1;n<t;++n)i>r.entry(s,n)&&(e[s]=n,i=r.entry(s,e[s]))}for(let s=0;s<t;++s)r.entry(s,e[s])<r.entry(a,e[a])&&(a=s);let h=e[a],l=s[a][0],_=s[h][0],c=l.isLeaf?[l.index]:l.index,u=_.isLeaf?[_.index]:_.index,d=c.concat(u),m=new Cluster(this._id++,l,_,r.entry(a,h),null,d);l.parent=m,_.parent=m,s[a].unshift(m),i[a]+=i[h];for(let e=0;e<t;++e){const t=r.entry(a,e),s=r.entry(h,e);let o;switch(n){case"single":o=Math.min(t,s);break;case"complete":o=Math.max(t,s);break;case"average":o=(i[a]*t+i[h]*s)/(i[a]+i[e])}r.set_entry(e,a,o),r.set_entry(a,e,o)}r.set_entry(a,a,1/0);for(let e=0;e<t;++e)r.set_entry(e,h,1/0),r.set_entry(h,e,1/0);
/* for (let j = 0; j < n; ++j) {
                if (d_min[j] === c2) {
                    d_min[j] = c1;
                }
                if (D.entry(c1, j) < D.entry(c1, d_min[c1])) {
                    d_min[c1] = j;
                }
            } */o=m}return o}}class Cluster{constructor(t,e,r,s,i,n,o,a){return this.id=t,this.left=e,this.right=r,this.dist=s,this.index=n,this.size=o??e.size+r.size,this.depth=a??1+Math.max(e.depth,r.depth),this.centroid=i??this._calculate_centroid(e,r),this.parent=null,this}_calculate_centroid(t,e){const r=t.size,s=e.size,i=t.centroid,n=e.centroid,o=this.size,a=t.centroid.length,h=new Float64Array(a);for(let t=0;t<a;++t)h[t]=(r*i[t]+s*n[t])/o;return h}get isLeaf(){return 0===this.depth}leaves(){if(this.isLeaf)return[this];const t=this.left,e=this.right;return(t.isLeaf?[t]:t.leaves()).concat(e.isLeaf?[e]:e.leaves())}descendants(){if(this.isLeaf)return[this];const t=this.left.descendants(),e=this.right.descendants();return t.concat(e).concat([this])}}
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
h.map(((t,e)=>[t,e])).filter((([t,e])=>t<i[e])).forEach((([t,e])=>{t<i[e]&&(i[e]=t,n[e]=a)}))}})),min(i)>=0)return!0;
// execute all improvements
for(;min(i)<0;){
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
for(let e=1;e<t;++e){let t=1/0;l=i.choice(s.filter((t=>a.findIndex((e=>e===t))<0)),n);for(let e=0;e<n;++e){let s=0;const i=l[e],o=r[i];for(let t=0;t<n;++t){if(t===e)continue;const n=l[t],h=r[n];let _=this._get_distance(i,n,o,h)-min(a.map((t=>this._get_distance(n,t,h))));_<0&&(s+=_)}
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
 * @extends DR
 */class LSP extends DR{
/**
     * Least Squares Projection.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias LSP
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.neighbors = Math.max(Math.floor(N / 10), 2)] - number of neighbors to consider.
     * @param {Number} [parameters.control_points = Math.ceil(Math.sqrt(N))] - number of controlpoints
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {LSP}
     * @see {@link https://ieeexplore.ieee.org/document/4378370}
     * @todo accept precomputed distance matrix.
     */
constructor(t,e){return super(t,{neighbors:void 0,control_points:void 0,d:2,metric:euclidean,seed:1212},e),this.parameter("neighbors",Math.min(e.neighbors??Math.max(Math.floor(this._N/10),2),this._N-1)),this.parameter("control_points",Math.min(e.control_points??Math.ceil(Math.sqrt(this._N)),this._N-1)),this._is_initialized=!1,this}
/**
     *
     * @param {DR} DR - method used for position control points.
     * @param {Object} DR_parameters - Object containing parameters for the DR method which projects the control points
     * @returns {LSP}
     */init(t=MDS,e={},r=BallTree){if(this._is_initialized)return this;const s=this.X,i=this._N,n=this.parameter("neighbors"),o=this.parameter("d"),a=this.parameter("seed"),h=this.parameter("metric");e=Object.assign({d:o,metric:h,seed:a},e);const l=this.parameter("control_points"),_=new KMedoids(s,l,null,h).get_clusters().medoids,c=new Matrix(l,i,"zeros");_.forEach(((t,e)=>{c.set_entry(e,t,1)}));const u=new t(Matrix.from(_.map((t=>s.row(t)))),e).transform(),d=s.to2dArray,m=new r(d,h),f=new Matrix(i,i,"I"),p=-1/n;d.forEach(((t,e)=>{for(const{index:r}of m.search(t,n).iterate())e!==r&&f.set_entry(e,r,p)}));const w=f.concat(c,"vertical"),g=new Matrix(i,o,"zeros").concat(u,"vertical");return this._A=w,this._b=g,this._is_initialized=!0,this}
/**
     * Computes the projection.
     * @returns {Matrix} Returns the projection.
     */transform(){this.check_init();const t=this._A,e=this._b,r=t.transDot(t),s=t.transDot(e);return this.Y=Matrix.solve_CG(r,s,this._randomizer),this.projection}}
/**
 * @class
 * @alias TopoMap
 * @memberof module:dimensionality_reduction
 * @extends DR
 */class TopoMap extends DR{
/**
     * TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias TopoMap
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Function} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {TopoMap}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
constructor(t,e){return super(t,{metric:euclidean,seed:1212},e),[this._N,this._D]=this.X.shape,this._distance_matrix=new Matrix(this._N,this._N,0),this}
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
     */init(){const{metric:t}=this._parameters;return this.Y=new Matrix(this._N,2,0),this._Emst=this._make_minimum_spanning_tree(t),this._is_initialized=!0,this}
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
     */transform(){this._is_initialized||this.init();const t=this._Emst,e=this.Y.to2dArray,r=new DisjointSet(e.map(((t,e)=>(t.i=e,t))));for(const[s,i,n]of t){const t=r.find(e[s]),o=r.find(e[i]);t!==o&&(this.__align_components(t,o,n),r.union(t,o))}return this.projection}*generator(){this._is_initialized||this.init();const t=this._Emst,e=this.Y.to2dArray,r=new DisjointSet(e.map(((t,e)=>(t.i=e,t))));for(const[s,i,n]of t){const t=r.find(e[s]),o=r.find(e[i]);t!==o&&(this.__align_components(t,o,n),r.union(t,o),yield this.projection)}return this.projection}}
/**
 * @class
 * @alias SAMMON
 * @extends DR
 */class SAMMON extends DR{
/**
     * SAMMON's Mapping
     * @constructor
     * @memberof module:dimensionality_reduction
     * @alias SAMMON
     * @param {Matrix} X - the high-dimensional data.
     * @param {Object} parameters - Object containing parameterization of the DR method.
     * @param {Number} [parameters.d = 2] - the dimensionality of the projection.
     * @param {Function|"precomputed"} [parameters.metric = euclidean] - the metric which defines the distance between two points.
     * @param {"PCA"|"MDS"|"random"} [parameters.init = "random"] - Either "PCA" or "MDS", with which SAMMON initialiates the projection. With "random" a random matrix gets used as starting point.
     * @param {Object} [parameters.init_parameters] - Parameters for the {@link init}-DR method.
     * @param {Number} [parameters.seed = 1212] - the seed for the random number generator.
     * @returns {SAMMON}
     * @see {@link https://arxiv.org/pdf/2009.01512.pdf}
     */
constructor(t,e){return super(t,{magic:.1,d:2,metric:euclidean,seed:1212,init_DR:"random",init_parameters:{}},e),this}
/**
     * initializes the projection.
     * @private
     */init(){const t=this.X.shape[0],{d:e,metric:r,init_DR:s,init_parameters:i}=this._parameters;if("random"===s){const r=this._randomizer;this.Y=new Matrix(t,e,(()=>r.random))}else{if(!["PCA","MDS"].includes(s))throw new Error('init_DR needs to be either "random" or a DR method!');this.Y=Matrix.from("PCA"==s?PCA.transform(this.X,i):MDS.transform(this.X,i))}return this.distance_matrix="precomputed"==r?Matrix.from(this.X):distance_matrix(this.X,r),this}
/**
     * Transforms the inputdata {@link X} to dimenionality 2.
     * @param {Number} [max_iter=200] - Maximum number of iteration steps.
     * @returns {Matrix|Array} - The projection of {@link X}.
     */transform(t=200){this._is_initialized||this.init();for(let e=0;e<t;++e)this._step();return this.projection}
/**
     * Transforms the inputdata {@link X} to dimenionality 2.
     * @param {Number} [max_iter=200] - Maximum number of iteration steps.
     * @returns {Generator} - A generator yielding the intermediate steps of the projection of {@link X}.
     */*generator(t=200){this._is_initialized||this.init();for(let e=0;e<t;++e)this._step(),yield this.projection;return this.projection}_step(){const t=this.parameter("magic"),e=this.distance_matrix,r=this.X.shape[0],{d:s,metric:i}=this._parameters;let n=this.Y,o=new Matrix(r,s,0),a=new Float64Array(s);for(let h=0;h<r;++h){let l=new Float64Array(s),_=new Float64Array(s);const c=n.row(h);for(let t=0;t<r;++t){if(h===t)continue;const r=n.row(t),o=new Float64Array(s);for(let t=0;t<s;++t)o[t]=c[t]-r[t];const a=i(c,r),u=e.entry(h,t),d=u-a,m=Math.max(u*a,.01);for(let t=0;t<s;++t)l[t]+=o[t]*d/m,_[t]+=(d-Math.pow(o[t],2)*(1+d/a)/a)/m}for(let e=0;e<s;++e){const r=n.entry(h,e)+(t*l[e]/Math.abs(_[e])||0);o.set_entry(h,e,r),a[e]+=r}}for(let t=0;t<s;++t)a[t]/=r;for(let t=0;t<r;++t)for(let e=0;e<s;++e)n.set_entry(t,e,o.entry(t,e)-a[e]);return n}}class SQDMDS extends DR{
/**
     * SQuadMDS: a lean Stochastic Quartet MDS improving global structure preservation in neighbor embedding like t-SNE and UMAP.
     * @constructor
     * @memberof module:dimensionality_reduction
     * @param {Matrix|Number[][]} X
     * @param {Object} [parameters]
     * @param {Number} [parameters.d=2]
     * @param {Function} [parameters.metric = euclidean]
     * @param {Number} [parameters.decay_start = 0.1] - Percentage of iterations using exaggeration phase. If random init: it is recommended to start the decay later to give the time for the global config to adjust with big steps.
     * @param {Number} [parameters.decay_cte = 0.34] - Controls the decay of the learning parameter.
     * @param {Object} [parameters.init_DR]
     * @returns {SQDMDS}
     * @see {@link https://arxiv.org/pdf/2202.12087.pdf}
     */
constructor(t,e){return super(t,{d:2,metric:euclidean,seed:1212,decay_start:.1,decay_cte:.34,// 0.34
init_DR:{type:"random"}},e),this}
/**
     * @private
     */init(){const t=this._N,e=this.parameter("d");
// initialize helpers.
this._add=this.__add(e),this._sub_div=this.__sub_div(e),this._minus=this.__minus(e),this._mult=this.__mult(e),this._LR_init=Math.max(2,.005*t),this._LR=this._LR_init,this._offset=-Math.exp(-1/this.parameter("decay_cte")),this._momentums=new Matrix(t,e,0),this._grads=new Matrix(t,e,0),this._indices=linspace(0,t-1);
// initialize projection.
const r=this._randomizer;this.Y=new Matrix(t,e,(()=>r.random-.5));
// preparing metric for optimization.
const s=this.parameter("metric");"precomputed"===s?(this._HD_metric=function(t,e,r){return r.entry(t,e)},this._HD_metric_exaggeration=function(t,e,r){return Math.pow(r.entry(t,e),2)}):(this._HD_metric=function(t,e,r){return s(r.row(t),r.row(e))},this._HD_metric_exaggeration=s==euclidean?function(t,e,r){return euclidean_squared(r.row(t),r.row(e))}:function(t,e,r){return Math.pow(s(r.row(t),r.row(e)),2)})}
/**
     * Computes the projection.
     * @param {Number} [iterations=500] - Number of iterations.
     * @returns {Matrix|Number[][]} the projection.
     */transform(t=500){this.check_init(),this._decay_start=Math.round(this.parameter("decay_start")*t);for(let e=0;e<t;++e)this._step(e,t);return this.projection}
/**
     * Computes the projection.
     * @param {Number} [iterations=500] - number of iterations.
     * @yields {Matrix|Number[][]} the intermediate steps of the projection.
     */*generator(t=500){this.check_init(),this._decay_start=Math.round(this.parameter("decay_start")*t);for(let e=0;e<t;++e)this._step(e,t),yield this.projection;return this.projection}
/**
     * Performs an optimization step.
     * @private
     * @param {Number} i - Acutal iteration.
     * @param {Number} iterations - Number of iterations.
     */_step(t,e){const r=this._decay_start;if(t>r){const s=this.parameter("decay_cte"),i=this._offset,n=(t-r)/(e-r);this._LR=this._LR_init*(Math.exp(-n*n/s)+i),this._distance_exaggeration=!1}else this._distance_exaggeration=!0;this._nestrov_iteration(this._distance_exaggeration)}
/**
     * Creates quartets of non overlapping indices.
     * @private
     * @returns {Number[][]}
     */__quartets(){const t=this._N,e=t-t%4,r=this._randomizer.choice(this._indices,e),s=[];for(let t=0;t<e;t+=4)s.push(Uint32Array.of(r[t],r[t+1],r[t+2],r[t+3]));return s}
/**
     * Computes and applies gradients, and updates momentum.
     * @private
     * @param {Boolean} distance_exaggeration
     */_nestrov_iteration(t){const e=this._momentums.mult(.99,{inline:!0}),r=this._LR,s=this._fill_MDS_grads(this.Y.add(e),this._grads,t),[i,n]=e.shape;for(let t=0;t<i;++t){const i=s.row(t),o=norm(i);if(0==o)continue;const a=r/o,h=e.row(t);for(let t=0;t<n;++t)h[t]-=a*i[t]}// momentums -= (LR / norm) * grads
this.Y.add(e,{inline:!0})}
/**
     * Computes the gradients.
     * @param {Matrix} Y - The Projection.
     * @param {Matrix} grads - The gradients.
     * @param {Boolean} [exaggeration = false] - Whether or not to use early exaggeration.
     * @param {Boolean} [zero_grad = true] - Whether or not to reset the gradient in the beginning.
     * @returns {Matrix} the gradients.
     */_fill_MDS_grads(t,e,r=!1,s=!0){s&&
// compute new gradients
e.values.fill(0);const i=this._add,n=this.X;let o;o=1==r?this._HD_metric_exaggeration:this._HD_metric;const a=new Float64Array(6),h=this.__quartets();for(const[r,s,l,_]of h){
// compute quartet's HD distances.
a[0]=o(r,s,n),a[1]=o(r,l,n),a[2]=o(r,_,n),a[3]=o(s,l,n),a[4]=o(s,_,n),a[5]=o(l,_,n);const h=neumair_sum(a);if(h>0)for(let t=0;t<6;++t)a[t]/=h,a[t]+=1e-11;const[c,u,d,m]=this._compute_quartet_grads(t,[r,s,l,_],a);
// add is inline, row acces the matrix
i(e.row(r),c),i(e.row(s),u),i(e.row(l),d),i(e.row(_),m)}return e}
/**
     * Quartet gradients for a projection.
     * @private
     * @param {Matrix} Y - The acutal projection.
     * @param {Number[]} quartet - The indices of the quartet.
     * @param {Number[]} D_hd - The high-dimensional distances of the quartet.
     * @returns {Number[][]} the gradients for the quartet.
     */_compute_quartet_grads(t,e,[r,s,i,n,o,a]){const[h,l,_,c]=e.map((e=>t.row(e))),u=euclidean(h,l)+1e-12,d=euclidean(h,_)+1e-12,m=euclidean(h,c)+1e-12,f=euclidean(l,_)+1e-12,p=euclidean(l,c)+1e-12,w=euclidean(_,c)+1e-12,g=neumair_sum([u,d,m,f,p,w]),[y,M,x,A]=this._ABCD_grads(h,l,_,c,u,d,m,f,p,w,r,g),[b,v,z,D]=this._ABCD_grads(h,_,l,c,d,u,m,f,w,p,s,g),[j,N,E,k]=this._ABCD_grads(h,c,_,l,m,d,u,w,p,f,i,g),[R,S,q,B]=this._ABCD_grads(l,_,h,c,f,u,p,d,w,m,n,g),[L,X,F,P]=this._ABCD_grads(l,c,h,_,p,u,f,m,w,d,o,g),[$,T,Y,C]=this._ABCD_grads(_,c,h,l,w,d,f,m,p,u,a,g),K=this._add;
// LD distances, add a small number just in case
return[K(y,b,j,q,F,Y),K(M,z,k,R,L,C),K(x,v,E,S,P,$),K(A,D,N,B,X,T)]}
/**
     * Gradients for one element of the loss function's sum.
     * @private
     */_ABCD_grads(t,e,r,s,i,n,o,a,h,l,_,c){const u=i/c,d=(_-u)/c*2,m=this._minus,f=this._add,p=this._mult,w=this._sub_div;return[p(m(p(f(w(t,e,i),w(t,r,n),w(t,s,o)),u),w(t,e,i)),d),p(m(p(f(w(e,t,i),w(e,r,a),w(e,s,h)),u),w(e,t,i)),d),p(f(w(r,t,n),w(r,e,a),w(r,s,l)),u*d),p(f(w(s,t,o),w(s,e,h),w(s,r,l)),u*d)]}
/**
     * Inline!
     */__minus(t){return(e,r)=>{for(let s=0;s<t;++s)e[s]-=r[s];return e}}
/**
     * Inline!
     */__add(t){return(...e)=>{const r=e.length,s=e[0];for(let i=1;i<r;++i){const r=e[i];for(let e=0;e<t;++e)s[e]+=r[e]}return s}}
/**
     * Inline!
     */__mult(t){return(e,r)=>{for(let s=0;s<t;++s)e[s]*=r;return e}}
/**
     * Creates a new array <code>(x - y) / div</code>
     */__sub_div(t){return(e,r,s)=>Float64Array.from({length:t},((t,i)=>(e[i]-r[i])/s))}}var t="0.6.3";export{BallTree,DisjointSet,FASTMAP,Heap,Hierarchical_Clustering,ISOMAP,KMeans,KMedoids,KNN,LDA,LLE,LSP,LTSA,MDS,Matrix,OPTICS,PCA,Randomizer,SAMMON,SQDMDS,TSNE,TopoMap,TriMap,UMAP,canberra,chebyshev,cosine,distance_matrix,euclidean,euclidean_squared,hamming,inner_product,jaccard,k_nearest_neighbors,kahan_sum,linspace,manhattan,max,min,neumair_sum,norm,normalize,powell,qr_gramschmidt as qr,qr_householder,simultaneous_poweriteration,sokal_michener,t as version,yule};
//# sourceMappingURL=druid.esm.js.map
