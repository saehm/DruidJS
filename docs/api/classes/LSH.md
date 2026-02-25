[@saehrimnir/druidjs](../globals.md) / LSH

# Class: LSH\<T\>

Defined in: knn/LSH.js:34

Locality Sensitive Hashing (LSH) for approximate nearest neighbor search.

LSH uses hash functions that map similar items to the same buckets with high probability.
This implementation uses Random Projection hashing (SimHash-style) which works well for
cosine similarity and Euclidean distance.

Key concepts:

- Multiple hash tables increase recall probability
- Each hash function projects data onto random hyperplanes
- Points on the same side of hyperplanes are hashed together
- Combines results from all tables for better accuracy

Best suited for:

- High-dimensional data where exact methods fail
- Approximate nearest neighbor needs
- Large datasets where linear scan is too slow
- When some false positives/negatives are acceptable

## Template

## See

[https://en.wikipedia.org/wiki/Locality-sensitive_hashing](https://en.wikipedia.org/wiki/Locality-sensitive_hashing)

## Extends

- `KNN`

## Type Parameters

### T

`T` _extends_ `number`[] \| `Float64Array`

## Constructors

### Constructor

> **new LSH**\<`T`\>(`elements`, `parameters?`): `LSH`\<`T`\>

Defined in: knn/LSH.js:41

Creates a new LSH index.

#### Parameters

##### elements

`T`[]

Elements to index

##### parameters?

[`ParametersLSH`](../interfaces/ParametersLSH.md) = `...`

Configuration parameters

#### Returns

`LSH`\<`T`\>

#### Overrides

`KNN.constructor`

## Properties

### \_dim

> **\_dim**: `number`

Defined in: knn/LSH.js:76

---

### \_elements

> **\_elements**: `T`[]

Defined in: [knn/KNN.js:14](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L14)

#### Inherited from

`KNN._elements`

---

### \_hashTables

> **\_hashTables**: `Map`\<`string`, `number`[]\>[]

Defined in: knn/LSH.js:64

---

### \_metric

> **\_metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: knn/LSH.js:56

---

### \_numHashFunctions

> **\_numHashFunctions**: `number`

Defined in: knn/LSH.js:58

---

### \_numHashTables

> **\_numHashTables**: `number`

Defined in: knn/LSH.js:57

---

### \_offsets

> **\_offsets**: `number`[][]

Defined in: knn/LSH.js:72

---

### \_parameters

> **\_parameters**: [`ParametersLSH`](../interfaces/ParametersLSH.md)

Defined in: [knn/KNN.js:16](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L16)

#### Inherited from

`KNN._parameters`

---

### \_projections

> **\_projections**: `Float64Array`\<`ArrayBufferLike`\>[][]

Defined in: knn/LSH.js:68

---

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: knn/LSH.js:60

---

### \_seed

> **\_seed**: `number`

Defined in: knn/LSH.js:59

---

### \_type

> **\_type**: `"array"` \| `"typed"`

Defined in: [knn/KNN.js:18](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L18)

#### Inherited from

`KNN._type`

## Methods

### add()

> **add**(`elements`): `LSH`\<`T`\>

Defined in: knn/LSH.js:169

Add elements to the LSH index.

#### Parameters

##### elements

`T`[]

#### Returns

`LSH`\<`T`\>

---

### search()

> **search**(`query`, `k?`): `object`[]

Defined in: knn/LSH.js:202

Search for k approximate nearest neighbors.

#### Parameters

##### query

`T`

##### k?

`number` = `5`

#### Returns

`object`[]

#### Overrides

`KNN.search`

---

### search_by_index()

> **search_by_index**(`i`, `k?`): `object`[]

Defined in: knn/LSH.js:288

#### Parameters

##### i

`number`

##### k?

`number` = `5`

#### Returns

`object`[]

#### Overrides

`KNN.search_by_index`
