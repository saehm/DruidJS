[@saehrimnir/druidjs](../globals.md) / NaiveKNN

# Class: NaiveKNN\<T\>

Defined in: knn/NaiveKNN.js:20

Naive KNN implementation using a distance matrix.

This implementation pre-computes the entire distance matrix and performs
an exhaustive search. Best suited for small datasets or when a distance
matrix is already available.

## Template

## Extends

- `KNN`

## Type Parameters

### T

`T` _extends_ `number`[] \| `Float64Array`

## Constructors

### Constructor

> **new NaiveKNN**\<`T`\>(`elements`, `parameters?`): `NaiveKNN`\<`T`\>

Defined in: knn/NaiveKNN.js:27

Generates a KNN list with given `elements`.

#### Parameters

##### elements

`T`[]

Elements which should be added to the KNN list

##### parameters?

[`ParametersNaiveKNN`](../interfaces/ParametersNaiveKNN.md) = `{}`

#### Returns

`NaiveKNN`\<`T`\>

#### Overrides

`KNN.constructor`

## Properties

### \_D

> **\_D**: [`Matrix`](Matrix.md)

Defined in: knn/NaiveKNN.js:33

---

### \_elements

> **\_elements**: `T`[]

Defined in: [knn/KNN.js:14](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L14)

#### Inherited from

`KNN._elements`

---

### \_parameters

> **\_parameters**: [`ParametersNaiveKNN`](../interfaces/ParametersNaiveKNN.md)

Defined in: [knn/KNN.js:16](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L16)

#### Inherited from

`KNN._parameters`

---

### \_type

> **\_type**: `"array"` \| `"typed"`

Defined in: [knn/KNN.js:18](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L18)

#### Inherited from

`KNN._type`

---

### KNN

> **KNN**: [`Heap`](Heap.md)\<\{ `index`: `number`; `value`: `number`; \}\>[]

Defined in: knn/NaiveKNN.js:42

## Methods

### search()

> **search**(`t`, `k?`): `object`[]

Defined in: knn/NaiveKNN.js:98

#### Parameters

##### t

`T`

Query element.

##### k?

`number` = `5`

Number of nearest neighbors to return. Default is `5`

#### Returns

`object`[]

- List consists of the `k` nearest neighbors.

#### Overrides

`KNN.search`

---

### search_by_index()

> **search_by_index**(`i`, `k?`): `object`[]

Defined in: knn/NaiveKNN.js:61

#### Parameters

##### i

`number`

##### k?

`number` = `5`

#### Returns

`object`[]

#### Overrides

`KNN.search_by_index`
