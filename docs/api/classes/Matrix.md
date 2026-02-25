[@saehrimnir/druidjs](../globals.md) / Matrix

# Class: Matrix

Defined in: [matrix/Matrix.js:11](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L11)

## Constructors

### Constructor

> **new Matrix**(`rows`, `cols`, `value?`): `Matrix`

Defined in: [matrix/Matrix.js:31](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L31)

Creates a new Matrix. Entries are stored in a Float64Array.

#### Parameters

##### rows

`number`

The amount of rows of the matrix.

##### cols

`number`

The amount of columns of the matrix.

##### value?

Can be a function with row and col as parameters, a number, or
  "zeros", "identity" or "I", or "center".

  - **function**: for each entry the function gets called with the parameters for the actual row and column.
  - **string**: allowed are

      - "zero", creates a zero matrix.
      - "identity" or "I", creates an identity matrix.
      - "center", creates an center matrix.
  - **number**: create a matrix filled with the given value.

`string` | `number` | `Accessor`

#### Returns

`Matrix`

#### Example

```ts
let A = new Matrix(10, 10, () => Math.random()); //creates a 10 times 10 random matrix. let B = new
Matrix(3, 3, "I"); // creates a 3 times 3 identity matrix.
```

## Properties

### \_cols

> **\_cols**: `number`

Defined in: [matrix/Matrix.js:33](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L33)

***

### \_data

> **\_data**: `Float64Array`\<`ArrayBufferLike`\>

Defined in: [matrix/Matrix.js:34](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L34)

***

### \_rows

> **\_rows**: `number`

Defined in: [matrix/Matrix.js:32](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L32)

## Accessors

### shape

#### Get Signature

> **get** **shape**(): `number`[]

Defined in: [matrix/Matrix.js:879](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L879)

Returns the number of rows and columns of the Matrix.

##### Returns

`number`[]

An Array in the form [rows, columns].

#### Set Signature

> **set** **shape**(`parameter`): `void`

Defined in: [matrix/Matrix.js:891](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L891)

Returns the matrix in the given shape with the given function which returns values for the entries of the matrix.

##### Parameters

###### parameter

\[`number`, `number`, `Accessor`\]

Takes an Array in the form [rows, cols, value], where rows and
  cols are the number of rows and columns of the matrix, and value is a function which takes two parameters (row
  and col) which has to return a value for the colth entry of the rowth row.

##### Returns

`void`

***

### T

#### Get Signature

> **get** **T**(): `Matrix`

Defined in: [matrix/Matrix.js:309](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L309)

Returns a new transposed Matrix. Short-form of `transpose`.

##### Returns

`Matrix`

***

### values

#### Get Signature

> **get** **values**(): `Float64Array`\<`ArrayBufferLike`\>

Defined in: [matrix/Matrix.js:970](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L970)

Returns the entries of the Matrix.

##### Returns

`Float64Array`\<`ArrayBufferLike`\>

## Methods

### \_apply()

> **\_apply**(`value`, `f`): `Matrix`

Defined in: [matrix/Matrix.js:725](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L725)

#### Parameters

##### value

`number` | `number`[] | `Float64Array`\<`ArrayBufferLike`\> | `Matrix`

##### f

(`d`, `v`) => `number`

#### Returns

`Matrix`

***

### \_apply\_colwise\_array()

> **\_apply\_colwise\_array**(`values`, `f`): `Matrix`

Defined in: [matrix/Matrix.js:708](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L708)

#### Parameters

##### values

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

##### f

(`d`, `v`) => `number`

#### Returns

`Matrix`

***

### \_apply\_rowwise\_array()

> **\_apply\_rowwise\_array**(`values`, `f`): `Matrix`

Defined in: [matrix/Matrix.js:700](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L700)

#### Parameters

##### values

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

##### f

(`d`, `v`) => `number`

#### Returns

`Matrix`

***

### \[iterator\]()

> **\[iterator\]**(): `Generator`\<`Float64Array`\<`ArrayBufferLike`\>, `void`, `unknown`\>

Defined in: [matrix/Matrix.js:178](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L178)

Makes a `Matrix` object an iterable object.

#### Returns

`Generator`\<`Float64Array`\<`ArrayBufferLike`\>, `void`, `unknown`\>

#### Yields

***

### add()

> **add**(`value`, `options?`): `Matrix`

Defined in: [matrix/Matrix.js:850](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L850)

Entrywise addition with `value`.

#### Parameters

##### value

`number` | `number`[] | `Float64Array`\<`ArrayBufferLike`\> | `Matrix`

##### options?

###### inline?

`boolean` = `false`

If true, applies addition to the element, otherwise it creates first a
  copy and applies the addition on the copy. Default is `false`

#### Returns

`Matrix`

#### Example

```ts
let A = Matrix.from([ [1, 2], [3, 4], ]); // a 2 by 2 matrix. let B = A.clone(); // B == A;

    A.add(2); // [[3, 4], [5, 6]];
    A.add(B); // [[2, 4], [6, 8]];
```

***

### add\_entry()

> **add\_entry**(`row`, `col`, `value`): `Matrix`

Defined in: [matrix/Matrix.js:275](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L275)

Adds a given [value](#add-entry) to the [col](#add-entry)<sup>th</sup> entry from the [row](#add-entry)<sup>th</sup> row of the
Matrix.

#### Parameters

##### row

`number`

##### col

`number`

##### value

`number`

#### Returns

`Matrix`

***

### asArray()

> **asArray**(): `number`[][]

Defined in: [matrix/Matrix.js:920](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L920)

Returns the Matrix as a Array of Arrays.

#### Returns

`number`[][]

***

### clone()

> **clone**(): `Matrix`

Defined in: [matrix/Matrix.js:788](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L788)

Clones the Matrix.

#### Returns

`Matrix`

***

### col()

> **col**(`col`): `Float64Array`\<`ArrayBufferLike`\>

Defined in: [matrix/Matrix.js:233](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L233)

Returns the col<sup>th</sup> column from the Matrix.

#### Parameters

##### col

`number`

#### Returns

`Float64Array`\<`ArrayBufferLike`\>

***

### concat()

> **concat**(`B`, `type?`): `Matrix`

Defined in: [matrix/Matrix.js:561](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L561)

Appends matrix `B` to the matrix.

#### Parameters

##### B

`Matrix`

Matrix to append.

##### type?

Type of concatenation. Default is
  `"horizontal"`

`"horizontal"` | `"vertical"` | `"diag"`

#### Returns

`Matrix`

#### Example

```ts
let A = Matrix.from([ [1, 1], [1, 1], ]); // 2 by 2 matrix filled with ones. let B = Matrix.from([ [2,
2], [2, 2], ]); // 2 by 2 matrix filled with twos.

    A.concat(B, "horizontal"); // 2 by 4 matrix. [[1, 1, 2, 2], [1, 1, 2, 2]]
    A.concat(B, "vertical"); // 4 by 2 matrix. [[1, 1], [1, 1], [2, 2], [2, 2]]
    A.concat(B, "diag"); // 4 by 4 matrix. [[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 2, 2], [0, 0, 2, 2]]
```

***

### diag()

> **diag**(): `Float64Array`\<`ArrayBufferLike`\>

Defined in: [matrix/Matrix.js:933](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L933)

Returns the diagonal of the Matrix.

#### Returns

`Float64Array`\<`ArrayBufferLike`\>

***

### divide()

> **divide**(`value`, `options?`): `Matrix`

Defined in: [matrix/Matrix.js:831](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L831)

Entrywise division with `value`.

#### Parameters

##### value

`number` | `number`[] | `Float64Array`\<`ArrayBufferLike`\> | `Matrix`

##### options?

###### inline?

`boolean` = `false`

If true, applies division to the element, otherwise it creates first a
  copy and applies the division on the copy. Default is `false`

#### Returns

`Matrix`

#### Example

```ts
let A = Matrix.from([ [1, 2], [3, 4], ]); // a 2 by 2 matrix. let B = A.clone(); // B == A;

    A.divide(2); // [[0.5, 1], [1.5, 2]];
    A.divide(B); // [[1, 1], [1, 1]];
```

***

### dot()

> **dot**(`B`): `Matrix`

Defined in: [matrix/Matrix.js:386](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L386)

Returns the dot product. If `B` is an Array or Float64Array then an Array gets returned. If `B` is a Matrix then
a Matrix gets returned.

#### Parameters

##### B

The right side

`number`[] | `Float64Array`\<`ArrayBufferLike`\> | `Matrix`

#### Returns

`Matrix`

***

### dotTrans()

> **dotTrans**(`B`): `Matrix`

Defined in: [matrix/Matrix.js:487](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L487)

Returns the dot product with the transposed version of `B`. If `B` is an Array or Float64Array then an Array gets
returned. If `B` is a Matrix then a Matrix gets returned.

#### Parameters

##### B

The right side

`number`[] | `Float64Array`\<`ArrayBufferLike`\> | `Matrix`

#### Returns

`Matrix`

***

### entry()

> **entry**(`row`, `col`): `number`

Defined in: [matrix/Matrix.js:248](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L248)

Returns the `col`<sup>th</sup> entry from the `row`<sup>th</sup> row of the Matrix.

#### Parameters

##### row

`number`

##### col

`number`

#### Returns

`number`

***

### gather()

> **gather**(`row_indices`, `col_indices`): `Matrix`

Defined in: [matrix/Matrix.js:658](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L658)

Returns a new array gathering entries defined by the indices given by argument.

#### Parameters

##### row\_indices

`number`[]

Array consists of indices of rows for gathering entries of this matrix

##### col\_indices

`number`[]

Array consists of indices of cols for gathering entries of this matrix

#### Returns

`Matrix`

***

### get\_block()

> **get\_block**(`start_row`, `start_col`, `end_row?`, `end_col?`): `Matrix`

Defined in: [matrix/Matrix.js:632](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L632)

Extracts the entries from the `start_row`<sup>th</sup> row to the `end_row`<sup>th</sup> row, the
`start_col`<sup>th</sup> column to the `end_col`<sup>th</sup> column of the matrix. If `end_row` or `end_col` is
empty, the respective value is set to `this.rows` or `this.cols`.

#### Parameters

##### start\_row

`number`

##### start\_col

`number`

##### end\_row?

`number` | `null`

##### end\_col?

`number` | `null`

#### Returns

`Matrix`

Returns a `end_row` - `start_row` times `end_col` - `start_col` matrix, with respective entries
  from the matrix.

#### Example

```ts
let A = Matrix.from([ [1, 2, 3], [4, 5, 6], [7, 8, 9], ]); // a 3 by 3 matrix.

    A.get_block(1, 1); // [[5, 6], [8, 9]]
    A.get_block(0, 0, 1, 1); // [[1]]
    A.get_block(1, 1, 2, 2); // [[5]]
    A.get_block(0, 0, 2, 2); // [[1, 2], [4, 5]]
```

***

### inverse()

> **inverse**(): `Matrix`

Defined in: [matrix/Matrix.js:318](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L318)

Returns the inverse of the Matrix.

#### Returns

`Matrix`

***

### iterate\_rows()

> **iterate\_rows**(): `Generator`\<`Float64Array`\<`ArrayBufferLike`\>, `void`, `unknown`\>

Defined in: [matrix/Matrix.js:164](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L164)

Returns an generator yielding each row of the Matrix.

#### Returns

`Generator`\<`Float64Array`\<`ArrayBufferLike`\>, `void`, `unknown`\>

#### Yields

***

### mean()

> **mean**(): `number`

Defined in: [matrix/Matrix.js:949](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L949)

Returns the mean of all entries of the Matrix.

#### Returns

`number`

***

### meanCols()

> **meanCols**(): `Float64Array`\<`ArrayBufferLike`\>

Defined in: [matrix/Matrix.js:1000](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L1000)

Returns the mean of each column of the matrix.

#### Returns

`Float64Array`\<`ArrayBufferLike`\>

***

### meanRows()

> **meanRows**(): `Float64Array`\<`ArrayBufferLike`\>

Defined in: [matrix/Matrix.js:980](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L980)

Returns the mean of each row of the matrix.

#### Returns

`Float64Array`\<`ArrayBufferLike`\>

***

### mult()

> **mult**(`value`, `options?`): `Matrix`

Defined in: [matrix/Matrix.js:812](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L812)

Entrywise multiplication with `value`.

#### Parameters

##### value

`number` | `number`[] | `Float64Array`\<`ArrayBufferLike`\> | `Matrix`

##### options?

###### inline?

`boolean` = `false`

If true, applies multiplication to the element, otherwise it creates
  first a copy and applies the multiplication on the copy. Default is `false`

#### Returns

`Matrix`

#### Example

```ts
let A = Matrix.from([ [1, 2], [3, 4], ]); // a 2 by 2 matrix. let B = A.clone(); // B == A;

    A.mult(2); // [[2, 4], [6, 8]];
    A.mult(B); // [[1, 4], [9, 16]];
```

***

### outer()

> **outer**(`B`): `Matrix`

Defined in: [matrix/Matrix.js:527](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L527)

Computes the outer product from `this` and `B`.

#### Parameters

##### B

`Matrix`

#### Returns

`Matrix`

***

### row()

> **row**(`row`): `Float64Array`\<`ArrayBufferLike`\>

Defined in: [matrix/Matrix.js:153](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L153)

Returns the `row`<sup>th</sup> row from the Matrix.

#### Parameters

##### row

`number`

#### Returns

`Float64Array`\<`ArrayBufferLike`\>

***

### set\_block()

> **set\_block**(`offset_row`, `offset_col`, `B`): `Matrix`

Defined in: [matrix/Matrix.js:602](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L602)

Writes the entries of B in A at an offset position given by `offset_row` and `offset_col`.

#### Parameters

##### offset\_row

`number`

##### offset\_col

`number`

##### B

`Matrix`

#### Returns

`Matrix`

***

### set\_entry()

> **set\_entry**(`row`, `col`, `value`): `Matrix`

Defined in: [matrix/Matrix.js:261](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L261)

Sets the [col](#set-entry)<sup>th</sup> entry from the [row](#set-entry)<sup>th</sup> row of the Matrix to the given
[value](#set-entry).

#### Parameters

##### row

`number`

##### col

`number`

##### value

`number`

#### Returns

`Matrix`

***

### set\_row()

> **set\_row**(`row`, `values`): `Matrix`

Defined in: [matrix/Matrix.js:191](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L191)

Sets the entries of `row`<sup>th</sup> row from the Matrix to the entries from `values`.

#### Parameters

##### row

`number`

##### values

`number`[]

#### Returns

`Matrix`

***

### sub()

> **sub**(`value`, `options?`): `Matrix`

Defined in: [matrix/Matrix.js:869](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L869)

Entrywise subtraction with `value`.

#### Parameters

##### value

`number` | `number`[] | `Float64Array`\<`ArrayBufferLike`\> | `Matrix`

##### options?

###### inline?

`boolean` = `false`

If true, applies subtraction to the element, otherwise it creates first
  a copy and applies the subtraction on the copy. Default is `false`

#### Returns

`Matrix`

#### Example

```ts
let A = Matrix.from([ [1, 2], [3, 4], ]); // a 2 by 2 matrix. let B = A.clone(); // B == A;

    A.sub(2); // [[-1, 0], [1, 2]];
    A.sub(B); // [[0, 0], [0, 0]];
```

***

### sub\_entry()

> **sub\_entry**(`row`, `col`, `value`): `Matrix`

Defined in: [matrix/Matrix.js:289](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L289)

Subtracts a given [value](#sub-entry) from the [col](#sub-entry)<sup>th</sup> entry from the [row](#sub-entry)<sup>th</sup> row of the
Matrix.

#### Parameters

##### row

`number`

##### col

`number`

##### value

`number`

#### Returns

`Matrix`

***

### sum()

> **sum**(): `number`

Defined in: [matrix/Matrix.js:960](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L960)

Returns the sum oof all entries of the Matrix.

#### Returns

`number`

***

### swap\_rows()

> **swap\_rows**(`row1`, `row2`): `Matrix`

Defined in: [matrix/Matrix.js:216](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L216)

Swaps the rows `row1` and `row2` of the Matrix.

#### Parameters

##### row1

`number`

##### row2

`number`

#### Returns

`Matrix`

***

### to2dArray()

> **to2dArray**(): `Float64Array`\<`ArrayBufferLike`\>[]

Defined in: [matrix/Matrix.js:907](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L907)

Returns the Matrix as a Array of Float64Arrays.

#### Returns

`Float64Array`\<`ArrayBufferLike`\>[]

***

### transDot()

> **transDot**(`B`): `Matrix`

Defined in: [matrix/Matrix.js:436](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L436)

Transposes the current matrix and returns the dot product with `B`. If `B` is an Array or Float64Array then an
Array gets returned. If `B` is a Matrix then a Matrix gets returned.

#### Parameters

##### B

The right side

`number`[] | `Float64Array`\<`ArrayBufferLike`\> | `Matrix`

#### Returns

`Matrix`

***

### transpose()

> **transpose**(): `Matrix`

Defined in: [matrix/Matrix.js:299](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L299)

Returns a new transposed Matrix.

#### Returns

`Matrix`

***

### det()

> `static` **det**(`A`): `number`

Defined in: [matrix/Matrix.js:1124](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L1124)

Computes the determinante of `A`, by using the `LU` decomposition of `A`.

#### Parameters

##### A

`Matrix`

#### Returns

`number`

The determinate of the Matrix `A`.

***

### from()

> `static` **from**(`A`): `Matrix`

Defined in: [matrix/Matrix.js:99](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L99)

Creates a Matrix out of `A`.

#### Parameters

##### A

The matrix, array, or number, which should converted to a Matrix.

`Matrix` | `Float64Array`\<`ArrayBufferLike`\>[] | `number`[][]

#### Returns

`Matrix`

#### Example

```ts
let A = Matrix.from([ [1, 0], [0, 1], ]); //creates a two by two identity matrix.
```

***

### from\_diag()

> `static` **from\_diag**(`v`): `Matrix`

Defined in: [matrix/Matrix.js:124](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L124)

Creates a Matrix with the diagonal being the values of `v`.

#### Parameters

##### v

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

#### Returns

`Matrix`

#### Example

```ts
let S = Matrix.from_diag([1, 2, 3]); // creates [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
```

***

### from\_vector()

> `static` **from\_vector**(`v`, `type`): `Matrix`

Defined in: [matrix/Matrix.js:138](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L138)

Creates a Matrix with the diagonal being the values of `v`.

#### Parameters

##### v

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

##### type

`"col"` | `"row"`

#### Returns

`Matrix`

#### Example

```ts
let S = Matrix.from_diag([1, 2, 3]); // creates [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
```

***

### is2dArray()

> `static` **is2dArray**(`A`): A is Float64Array\<ArrayBufferLike\>\[\] \| number\[\]\[\]

Defined in: [matrix/Matrix.js:1190](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L1190)

#### Parameters

##### A

`any`[]

#### Returns

A is Float64Array\<ArrayBufferLike\>\[\] \| number\[\]\[\]

***

### isArray()

> `static` **isArray**(`A`): A is unknown\[\] \| number\[\] \| Float32Array\<ArrayBufferLike\> \| Float64Array\<ArrayBufferLike\>

Defined in: [matrix/Matrix.js:1182](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L1182)

#### Parameters

##### A

`unknown`

#### Returns

A is unknown\[\] \| number\[\] \| Float32Array\<ArrayBufferLike\> \| Float64Array\<ArrayBufferLike\>

***

### LU()

> `static` **LU**(`A`): `object`

Defined in: [matrix/Matrix.js:1090](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L1090)

`LU` decomposition of the Matrix `A`. Creates two matrices, so that the dot product `LU` equals `A`.

#### Parameters

##### A

`Matrix`

#### Returns

`object`

The left triangle matrix `L` and the upper triangle matrix `U`.

##### L

> **L**: `Matrix`

##### U

> **U**: `Matrix`

***

### solve()

> `static` **solve**(`A`, `b`): `Matrix`

Defined in: [matrix/Matrix.js:1060](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L1060)

Solves the equation `Ax = b`. Returns the result `x`.

#### Parameters

##### A

Matrix or LU Decomposition

`Matrix` | \{ `L`: `Matrix`; `U`: `Matrix`; \}

##### b

`Matrix`

Matrix

#### Returns

`Matrix`

***

### solve\_CG()

> `static` **solve\_CG**(`A`, `b`, `randomizer?`, `tol?`): `Matrix`

Defined in: [matrix/Matrix.js:1024](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L1024)

Solves the equation `Ax = b` using the conjugate gradient method. Returns the result `x`.

#### Parameters

##### A

`Matrix`

Matrix

##### b

`Matrix`

Matrix

##### randomizer?

[`Randomizer`](Randomizer.md) | `null`

##### tol?

`number` = `1e-3`

Default is `1e-3`

#### Returns

`Matrix`

***

### SVD()

> `static` **SVD**(`M`, `k?`): `object`

Defined in: [matrix/Matrix.js:1160](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/matrix/Matrix.js#L1160)

Computes the `k` components of the SVD decomposition of the matrix `M`.

#### Parameters

##### M

`Matrix`

##### k?

`number` = `2`

Default is `2`

#### Returns

`object`

##### Sigma

> **Sigma**: `Float64Array`

##### U

> **U**: `Float64Array`\<`ArrayBufferLike`\>[]

##### V

> **V**: `Float64Array`\<`ArrayBufferLike`\>[]
