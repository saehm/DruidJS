
<a href="#module_matrix" name="module_matrix">#</a> <code>**matrix**</code>



* [matrix](#module_matrix)

    * [distance_matrix(A, [metric])](#distance_matrix)

    * [k_nearest_neigbhors(A, k, [metric])](#k_nearest_neigbhors)

    * [linspace(start, end, [number])](#linspace)

    * [norm(v, [metric])](#norm)

    * [normalize(v, metric)](#normalize)

    * _global_
        * [Matrix](#Matrix)

            * [new exports.Matrix(rows, cols, value)](#new_Matrix_new)

            * _instance_
                * [.T](#Matrix+T)

                * [.shape](#Matrix+shape)

                * [.shape](#Matrix+shape)

                * [.to2dArray](#Matrix+to2dArray)

                * [.asArray](#Matrix+asArray)

                * [.diag](#Matrix+diag)

                * [.mean](#Matrix+mean)

                * [.sum](#Matrix+sum)

                * [.values](#Matrix+values)

                * [.meanRows](#Matrix+meanRows)

                * [.meanCols](#Matrix+meanCols)

                * [.row(row)](#Matrix+row)

                * [.iterate_rows()](#Matrix+iterate_rows)

                * [.set_row(row, values)](#Matrix+set_row)

                * [.col(col)](#Matrix+col)

                * [.entry(row, col)](#Matrix+entry)

                * [.set_entry(row, col, value)](#Matrix+set_entry)

                * [.transpose()](#Matrix+transpose)

                * [.inverse()](#Matrix+inverse)

                * [.dot(B)](#Matrix+dot)

                * [.outer(B)](#Matrix+outer)

                * [.concat(B, [type])](#Matrix+concat)

                * [.set_block(offset_row, offset_col, B)](#Matrix+set_block)

                * [.get_block(start_row, start_col, [end_row], [end_col])](#Matrix+get_block)

                * [.gather(row_indices, col_indices)](#Matrix+gather)

                * [.clone()](#Matrix+clone)

                * [.mult(value)](#Matrix+mult)

                * [.divide(value)](#Matrix+divide)

                * [.add(value)](#Matrix+add)

                * [.sub(value)](#Matrix+sub)

            * _static_
                * [.from(A, [type])](#Matrix.from)

                * [.solve_CG(A, b, [randomizer], [tol])](#Matrix.solve_CG)

                * [.solve(A, b)](#Matrix.solve)

                * [.LU(A)](#Matrix.LU)

                * [.det(A)](#Matrix.det)

                * [.SVD(M, [k])](#Matrix.SVD)



<a href="#distance_matrix" name="distance_matrix">#</a> <code>*matrix***distance_matrix**</code>
(A, [metric])

Computes the distance matrix of datamatrix [A](A).


- A [<code>Matrix</code>](#Matrix) - Matrix.
- [metric] <code>function</code> <code> = euclidean</code> - The diistance metric.

**Returns**: [<code>Matrix</code>](#Matrix) - D - The distance matrix of [A](A).  

<a href="#k_nearest_neigbhors" name="k_nearest_neigbhors">#</a> <code>*matrix***k_nearest_neigbhors**</code>
(A, k, [metric])

Computes the k-nearest neighbors of each row of [A](A).


- A [<code>Matrix</code>](#Matrix) - Either the data matrix, or a distance matrix.
- k <code>Number</code> - The number of neighbors to compute.
- [metric] <code>function</code> | <code>&quot;precomputed&quot;</code> <code> = euclidean</code>

**Returns**: <code>Array.&lt;Object&gt;</code> - -  

<a href="#linspace" name="linspace">#</a> <code>*matrix***linspace**</code>
(start, end, [number])

Creates an Array containing [number](number) numbers from [start](start) to [end](end).
If <code>[number](number) = null</null>.


- start <code>Number</code> - Start value.
- end <code>Number</code> - End value.
- [number] <code>Number</code> <code> = </code> - Number of number between [start](start) and [end](end).

**Returns**: <code>Array</code> - - An array with [number](number) entries, beginning at [start](start) ending at [end](end).  

<a href="#norm" name="norm">#</a> <code>*matrix***norm**</code>
(v, [metric])

Computes the norm of a vector, by computing its distance to **0**.


- v [<code>Matrix</code>](#Matrix) | <code>Array.&lt;Number&gt;</code> | <code>Float64Array</code> - Vector.
- [metric] <code>function</code> <code> = euclidean</code> - Which metric should be used to compute the norm.

**Returns**: <code>Number</code> - - The norm of [v](v).  

<a href="#normalize" name="normalize">#</a> <code>*matrix***normalize**</code>
(v, metric)

Normalizes Vector [v](v).


- v <code>Array.&lt;Number&gt;</code> | <code>Float64Array</code> - Vector
- metric <code>function</code>

**Returns**: <code>Array.&lt;Number&gt;</code> \| <code>Float64Array</code> - - The normalized vector with length 1.  

<a href="#Matrix" name="Matrix">#</a> <code>*matrix***Matrix**</code>



* [Matrix](#Matrix)

    * [new exports.Matrix(rows, cols, value)](#new_Matrix_new)

    * _instance_
        * [.T](#Matrix+T)

        * [.shape](#Matrix+shape)

        * [.shape](#Matrix+shape)

        * [.to2dArray](#Matrix+to2dArray)

        * [.asArray](#Matrix+asArray)

        * [.diag](#Matrix+diag)

        * [.mean](#Matrix+mean)

        * [.sum](#Matrix+sum)

        * [.values](#Matrix+values)

        * [.meanRows](#Matrix+meanRows)

        * [.meanCols](#Matrix+meanCols)

        * [.row(row)](#Matrix+row)

        * [.iterate_rows()](#Matrix+iterate_rows)

        * [.set_row(row, values)](#Matrix+set_row)

        * [.col(col)](#Matrix+col)

        * [.entry(row, col)](#Matrix+entry)

        * [.set_entry(row, col, value)](#Matrix+set_entry)

        * [.transpose()](#Matrix+transpose)

        * [.inverse()](#Matrix+inverse)

        * [.dot(B)](#Matrix+dot)

        * [.outer(B)](#Matrix+outer)

        * [.concat(B, [type])](#Matrix+concat)

        * [.set_block(offset_row, offset_col, B)](#Matrix+set_block)

        * [.get_block(start_row, start_col, [end_row], [end_col])](#Matrix+get_block)

        * [.gather(row_indices, col_indices)](#Matrix+gather)

        * [.clone()](#Matrix+clone)

        * [.mult(value)](#Matrix+mult)

        * [.divide(value)](#Matrix+divide)

        * [.add(value)](#Matrix+add)

        * [.sub(value)](#Matrix+sub)

    * _static_
        * [.from(A, [type])](#Matrix.from)

        * [.solve_CG(A, b, [randomizer], [tol])](#Matrix.solve_CG)

        * [.solve(A, b)](#Matrix.solve)

        * [.LU(A)](#Matrix.LU)

        * [.det(A)](#Matrix.det)

        * [.SVD(M, [k])](#Matrix.SVD)



<a href="#new_Matrix_new" name="new_Matrix_new">#</a> new <code>**exports.Matrix**</code>
(rows, cols, value)

creates a new Matrix. Entries are stored in a Float64Array.


- rows <code>number</code> <code> = </code> - The amount of rows of the matrix.
- cols <code>number</code> <code> = </code> - The amount of columns of the matrix.
- value <code>function</code> | <code>string</code> | <code>number</code> <code> = 0</code> - Can be a function with row and col as parameters, a number, or "zeros", "identity" or "I", or "center".
 - **function**: for each entry the function gets called with the parameters for the actual row and column.
 - **string**: allowed are
     - "zero", creates a zero matrix.
     - "identity" or "I", creates an identity matrix.
     - "center", creates an center matrix.
 - **number**: create a matrix filled with the given value.

**Returns**: [<code>Matrix</code>](#Matrix) - returns a [rows](rows) times [cols](cols) Matrix filled with [value](value).  
**Example**  
```js
let A = new Matrix(10, 10, () => Math.random()); //creates a 10 times 10 random matrix.
let B = new Matrix(3, 3, "I"); // creates a 3 times 3 identity matrix.
```

<a href="#Matrix+T" name="Matrix+T">#</a> <code>*matrix*.**T**</code>


Returns a new transposed Matrix. Short-form of {@function transpose}.


<a href="#Matrix+shape" name="Matrix+shape">#</a> <code>*matrix*.**shape**</code>


Returns the number of rows and columns of the Matrix.

**Returns**: <code>Array</code> - An Array in the form [rows, columns].  

<a href="#Matrix+shape" name="Matrix+shape">#</a> <code>*matrix*.**shape**</code>


Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.


- parameter <code>Array</code> - takes an Array in the form [rows, cols, value], where rows and cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row and col) which has to return a value for the colth entry of the rowth row.


<a href="#Matrix+to2dArray" name="Matrix+to2dArray">#</a> <code>*matrix*.**to2dArray**</code>


Returns the Matrix as a Array of Float64Arrays.


<a href="#Matrix+asArray" name="Matrix+asArray">#</a> <code>*matrix*.**asArray**</code>


Returns the Matrix as a Array of Arrays.


<a href="#Matrix+diag" name="Matrix+diag">#</a> <code>*matrix*.**diag**</code>


Returns the diagonal of the Matrix.


<a href="#Matrix+mean" name="Matrix+mean">#</a> <code>*matrix*.**mean**</code>


Returns the mean of all entries of the Matrix.


<a href="#Matrix+sum" name="Matrix+sum">#</a> <code>*matrix*.**sum**</code>


Returns the sum oof all entries of the Matrix.


<a href="#Matrix+values" name="Matrix+values">#</a> <code>*matrix*.**values**</code>


Returns the sum oof all entries of the Matrix.


<a href="#Matrix+meanRows" name="Matrix+meanRows">#</a> <code>*matrix*.**meanRows**</code>


Returns the mean of each row of the matrix.


<a href="#Matrix+meanCols" name="Matrix+meanCols">#</a> <code>*matrix*.**meanCols**</code>


Returns the mean of each column of the matrix.


<a href="#Matrix+row" name="Matrix+row">#</a> <code>*matrix*.**row**</code>
(row)

Returns the [row](row)<sup>th</sup> row from the Matrix.


- row <code>Number</code>


<a href="#Matrix+iterate_rows" name="Matrix+iterate_rows">#</a> <code>*matrix*.**iterate_rows**</code>
()

Returns an generator yielding each row of the Matrix.


<a href="#Matrix+set_row" name="Matrix+set_row">#</a> <code>*matrix*.**set_row**</code>
(row, values)

Sets the entries of [row](row)<sup>th</sup> row from the Matrix to the entries from [values](values).


- row <code>int</code>
- values <code>Array</code>


<a href="#Matrix+col" name="Matrix+col">#</a> <code>*matrix*.**col**</code>
(col)

Returns the [col](col)<sup>th</sup> column from the Matrix.


- col <code>int</code>


<a href="#Matrix+entry" name="Matrix+entry">#</a> <code>*matrix*.**entry**</code>
(row, col)

Returns the [col](col)<sup>th</sup> entry from the [row](row)<sup>th</sup> row of the Matrix.


- row <code>int</code>
- col <code>int</code>


<a href="#Matrix+set_entry" name="Matrix+set_entry">#</a> <code>*matrix*.**set_entry**</code>
(row, col, value)

Sets the [col](col)<sup>th</sup> entry from the [row](row)<sup>th</sup> row of the Matrix to the given [value](value).


- row <code>int</code>
- col <code>int</code>
- value <code>float64</code>


<a href="#Matrix+transpose" name="Matrix+transpose">#</a> <code>*matrix*.**transpose**</code>
()

Returns a new transposed Matrix.


<a href="#Matrix+inverse" name="Matrix+inverse">#</a> <code>*matrix*.**inverse**</code>
()

Returns the inverse of the Matrix.


<a href="#Matrix+dot" name="Matrix+dot">#</a> <code>*matrix*.**dot**</code>
(B)

Returns the dot product. If [B](B) is an Array or Float64Array then an Array gets returned. If [B](B) is a Matrix then a Matrix gets returned.


- B [<code>Matrix</code>](#Matrix) | <code>Array</code> | <code>Float64Array</code> - the right side


<a href="#Matrix+outer" name="Matrix+outer">#</a> <code>*matrix*.**outer**</code>
(B)

Computes the outer product from [this](this) and [B](B).


- B [<code>Matrix</code>](#Matrix)


<a href="#Matrix+concat" name="Matrix+concat">#</a> <code>*matrix*.**concat**</code>
(B, [type])

Appends matrix [B](B) to the matrix.


- B [<code>Matrix</code>](#Matrix) - matrix to append.
- [type] <code>&quot;horizontal&quot;</code> | <code>&quot;vertical&quot;</code> | <code>&quot;diag&quot;</code> <code> = &quot;horizontal&quot;</code> - type of concatenation.

**Example**  
```js
let A = Matrix.from([[1, 1], [1, 1]]); // 2 by 2 matrix filled with ones.
let B = Matrix.from([[2, 2], [2, 2]]); // 2 by 2 matrix filled with twos.

A.concat(B, "horizontal"); // 2 by 4 matrix. [[1, 1, 2, 2], [1, 1, 2, 2]]
A.concat(B, "vertical"); // 4 by 2 matrix. [[1, 1], [1, 1], [2, 2], [2, 2]]
A.concat(B, "diag"); // 4 by 4 matrix. [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2], [0, 0, 2, 2]]
```

<a href="#Matrix+set_block" name="Matrix+set_block">#</a> <code>*matrix*.**set_block**</code>
(offset_row, offset_col, B)

Writes the entries of B in A at an offset position given by [offset_row](offset_row) and [offset_col](offset_col).


- offset_row <code>int</code>
- offset_col <code>int</code>
- B [<code>Matrix</code>](#Matrix)


<a href="#Matrix+get_block" name="Matrix+get_block">#</a> <code>*matrix*.**get_block**</code>
(start_row, start_col, [end_row], [end_col])

Extracts the entries from the [start_row](start_row)<sup>th</sup> row to the [end_row](end_row)<sup>th</sup> row, the [start_col](start_col)<sup>th</sup> column to the [end_col](end_col)<sup>th</sup> column of the matrix.
If [end_row](end_row) or [end_col](end_col) is empty, the respective value is set to [this.rows](this.rows) or [this.cols](this.cols).


- start_row <code>Number</code>
- start_col <code>Number</code>
- [end_row] <code>Number</code> <code> = </code>
- [end_col] <code>Number</code> <code> = </code>

**Returns**: [<code>Matrix</code>](#Matrix) - Returns a end_row - start_row times end_col - start_col matrix, with respective entries from the matrix.  
**Example**  
```js
let A = Matrix.from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]); // a 3 by 3 matrix.

A.get_block(1, 1); // [[5, 6], [8, 9]]
A.get_block(0, 0, 1, 1); // [[1]]
A.get_block(1, 1, 2, 2); // [[5]]
A.get_block(0, 0, 2, 2); // [[1, 2], [4, 5]]
```

<a href="#Matrix+gather" name="Matrix+gather">#</a> <code>*matrix*.**gather**</code>
(row_indices, col_indices)

Returns a new array gathering entries defined by the indices given by argument.


- row_indices <code>Array.&lt;Number&gt;</code> - Array consists of indices of rows for gathering entries of this matrix
- col_indices <code>Array.&lt;Number&gt;</code> - Array consists of indices of cols for gathering entries of this matrix


<a href="#Matrix+clone" name="Matrix+clone">#</a> <code>*matrix*.**clone**</code>
()

Clones the Matrix.


<a href="#Matrix+mult" name="Matrix+mult">#</a> <code>*matrix*.**mult**</code>
(value)

Entrywise multiplication with [value](value).


- value [<code>Matrix</code>](#Matrix) | <code>Array</code> | <code>Number</code>

**Example**  
```js
let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
let B = A.clone(); // B == A;

A.mult(2); // [[2, 4], [6, 8]];
A.mult(B); // [[1, 4], [9, 16]];
```

<a href="#Matrix+divide" name="Matrix+divide">#</a> <code>*matrix*.**divide**</code>
(value)

Entrywise division with [value](value).


- value [<code>Matrix</code>](#Matrix) | <code>Array</code> | <code>Number</code>

**Example**  
```js
let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
let B = A.clone(); // B == A;

A.divide(2); // [[0.5, 1], [1.5, 2]];
A.divide(B); // [[1, 1], [1, 1]];
```

<a href="#Matrix+add" name="Matrix+add">#</a> <code>*matrix*.**add**</code>
(value)

Entrywise addition with [value](value).


- value [<code>Matrix</code>](#Matrix) | <code>Array</code> | <code>Number</code>

**Example**  
```js
let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
let B = A.clone(); // B == A;

A.add(2); // [[3, 4], [5, 6]];
A.add(B); // [[2, 4], [6, 8]];
```

<a href="#Matrix+sub" name="Matrix+sub">#</a> <code>*matrix*.**sub**</code>
(value)

Entrywise subtraction with [value](value).


- value [<code>Matrix</code>](#Matrix) | <code>Array</code> | <code>Number</code>

**Example**  
```js
let A = Matrix.from([[1, 2], [3, 4]]); // a 2 by 2 matrix.
let B = A.clone(); // B == A;

A.sub(2); // [[-1, 0], [1, 2]];
A.sub(B); // [[0, 0], [0, 0]];
```

<a href="#Matrix.from" name="Matrix.from">#</a> <code>*Matrix*.**from**</code>
(A, [type])

Creates a Matrix out of [A](A).


- A [<code>Matrix</code>](#Matrix) | <code>Array</code> | <code>Float64Array</code> | <code>number</code> - The matrix, array, or number, which should converted to a Matrix.
- [type] <code>&quot;row&quot;</code> | <code>&quot;col&quot;</code> | <code>&quot;diag&quot;</code> <code> = &quot;row&quot;</code> - If [A](A) is a Array or Float64Array, then type defines if it is a row- or a column vector.

**Example**  
```js
let A = Matrix.from([[1, 0], [0, 1]]); //creates a two by two identity matrix.
let S = Matrix.from([1, 2, 3], "diag"); // creates a 3 by 3 matrix with 1, 2, 3 on its diagonal. [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
```

<a href="#Matrix.solve_CG" name="Matrix.solve_CG">#</a> <code>*Matrix*.**solve_CG**</code>
(A, b, [randomizer], [tol])

Solves the equation [A](A)x = [b](b) using the conjugate gradient method. Returns the result x.


- A [<code>Matrix</code>](#Matrix) - Matrix
- b [<code>Matrix</code>](#Matrix) - Matrix
- [randomizer] <code>Randomizer</code> <code> = </code>
- [tol] <code>Number</code> <code> = 1e-3</code>


<a href="#Matrix.solve" name="Matrix.solve">#</a> <code>*Matrix*.**solve**</code>
(A, b)

Solves the equation [A](A)x = [b](b). Returns the result x.


- A [<code>Matrix</code>](#Matrix) - Matrix or LU Decomposition
- b [<code>Matrix</code>](#Matrix) - Matrix


<a href="#Matrix.LU" name="Matrix.LU">#</a> <code>*Matrix*.**LU**</code>
(A)

[L](L)[U](U) decomposition of the Matrix [A](A). Creates two matrices, so that the dot product LU equals A.


- A [<code>Matrix</code>](#Matrix)

**Returns**: <code>Object</code> - result - Returns the left triangle matrix [L](L) and the upper triangle matrix [U](U).  

<a href="#Matrix.det" name="Matrix.det">#</a> <code>*Matrix*.**det**</code>
(A)

Computes the determinante of [A](A), by using the LU decomposition of [A](A).


- A [<code>Matrix</code>](#Matrix)

**Returns**: <code>Number</code> - det - Returns the determinate of the Matrix [A](A).  

<a href="#Matrix.SVD" name="Matrix.SVD">#</a> <code>*Matrix*.**SVD**</code>
(M, [k])

Computes the [k](k) components of the SVD decomposition of the matrix [M](M)


- M [<code>Matrix</code>](#Matrix)
- [k] <code>int</code> <code> = 2</code>

