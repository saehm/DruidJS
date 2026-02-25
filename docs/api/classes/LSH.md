[@saehrimnir/druidjs](../globals.md) / LSH

# Class: LSH\<T\>

Defined in: [knn/LSH.js:34](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L34)

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

[https://en.wikipedia.org/wiki/Locality-sensitive\_hashing](https://en.wikipedia.org/wiki/Locality-sensitive_hashing)

## Extends

- `KNN`

## Type Parameters

### T

`T` *extends* `number`[] \| `Float64Array`

## Constructors

### Constructor

> **new LSH**\<`T`\>(`elements`, `parameters?`): `LSH`\<`T`\>

Defined in: [knn/LSH.js:41](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L41)

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

Defined in: [knn/LSH.js:76](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L76)

***

### \_elements

> **\_elements**: `T`[]

Defined in: [knn/KNN.js:14](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L14)

#### Inherited from

`KNN._elements`

***

### \_hashTables

> **\_hashTables**: `Map`\<`string`, `number`[]\>[]

Defined in: [knn/LSH.js:64](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L64)

***

### \_metric

> **\_metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [knn/LSH.js:56](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L56)

***

### \_numHashFunctions

> **\_numHashFunctions**: `number`

Defined in: [knn/LSH.js:58](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L58)

***

### \_numHashTables

> **\_numHashTables**: `number`

Defined in: [knn/LSH.js:57](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L57)

***

### \_offsets

> **\_offsets**: `number`[][]

Defined in: [knn/LSH.js:72](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L72)

***

### \_parameters

> **\_parameters**: [`ParametersLSH`](../interfaces/ParametersLSH.md)

Defined in: [knn/KNN.js:16](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L16)

#### Inherited from

`KNN._parameters`

***

### \_projections

> **\_projections**: `Float64Array`\<`ArrayBufferLike`\>[][]

Defined in: [knn/LSH.js:68](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L68)

***

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: [knn/LSH.js:60](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L60)

***

### \_seed

> **\_seed**: `number`

Defined in: [knn/LSH.js:59](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L59)

***

### \_type

> **\_type**: `"array"` \| `"typed"`

Defined in: [knn/KNN.js:18](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L18)

#### Inherited from

`KNN._type`

## Methods

### add()

> **add**(`elements`): `LSH`\<`T`\>

Defined in: [knn/LSH.js:169](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L169)

Add elements to the LSH index.

#### Parameters

##### elements

`T`[]

#### Returns

`LSH`\<`T`\>

***

### search()

> **search**(`query`, `k?`): `object`[]

Defined in: [knn/LSH.js:202](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L202)

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

***

### search\_by\_index()

> **search\_by\_index**(`i`, `k?`): `object`[]

Defined in: [knn/LSH.js:288](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/LSH.js#L288)

#### Parameters

##### i

`number`

##### k?

`number` = `5`

#### Returns

`object`[]

#### Overrides

`KNN.search_by_index`
