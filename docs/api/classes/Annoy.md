[@saehrimnir/druidjs](../globals.md) / Annoy

# Class: Annoy\<T\>

Defined in: [knn/Annoy.js:45](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L45)

Annoy-style (Approximate Nearest Neighbors Oh Yeah) implementation using Random Projection Trees.

This implementation builds multiple random projection trees where each tree randomly selects
two points and splits the space based on a hyperplane equidistant between them.

Key features:
- Multiple random projection trees for better recall
- Each tree uses random hyperplanes for splitting
- Priority queue search for better recall
- Combines results from all trees

Best suited for:
- High-dimensional data
- Approximate nearest neighbor search
- Large datasets
- When high recall is needed with approximate methods

## Template

## See

 - [https://github.com/spotify/annoy](https://github.com/spotify/annoy)
 - [https://erikbern.com/2015/09/24/nearest-neighbors-and-vector-models-epilogue-curse-of-dimensionality.html](https://erikbern.com/2015/09/24/nearest-neighbors-and-vector-models-epilogue-curse-of-dimensionality.html)

## Extends

- `KNN`

## Type Parameters

### T

`T` *extends* `number`[] \| `Float64Array`

## Constructors

### Constructor

> **new Annoy**\<`T`\>(`elements`, `parameters?`): `Annoy`\<`T`\>

Defined in: [knn/Annoy.js:52](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L52)

Creates a new Annoy-style index with random projection trees.

#### Parameters

##### elements

`T`[]

Elements to index

##### parameters?

[`ParametersAnnoy`](../interfaces/ParametersAnnoy.md) = `...`

Configuration parameters

#### Returns

`Annoy`\<`T`\>

#### Overrides

`KNN.constructor`

## Properties

### \_elements

> **\_elements**: `T`[]

Defined in: [knn/KNN.js:14](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L14)

#### Inherited from

`KNN._elements`

***

### \_maxPointsPerLeaf

> **\_maxPointsPerLeaf**: `number`

Defined in: [knn/Annoy.js:69](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L69)

***

### \_metric

> **\_metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [knn/Annoy.js:67](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L67)

***

### \_numTrees

> **\_numTrees**: `number`

Defined in: [knn/Annoy.js:68](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L68)

***

### \_parameters

> **\_parameters**: [`ParametersAnnoy`](../interfaces/ParametersAnnoy.md)

Defined in: [knn/KNN.js:16](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L16)

#### Inherited from

`KNN._parameters`

***

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: [knn/Annoy.js:71](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L71)

***

### \_seed

> **\_seed**: `number`

Defined in: [knn/Annoy.js:70](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L70)

***

### \_type

> **\_type**: `"array"` \| `"typed"`

Defined in: [knn/KNN.js:18](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L18)

#### Inherited from

`KNN._type`

## Accessors

### num\_nodes

#### Get Signature

> **get** **num\_nodes**(): `number`

Defined in: [knn/Annoy.js:101](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L101)

Get the total number of nodes in all trees.

##### Returns

`number`

***

### num\_trees

#### Get Signature

> **get** **num\_trees**(): `number`

Defined in: [knn/Annoy.js:93](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L93)

Get the number of trees in the index.

##### Returns

`number`

## Methods

### add()

> **add**(`elements`): `Annoy`\<`T`\>

Defined in: [knn/Annoy.js:124](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L124)

Add elements to the Annoy index.

#### Parameters

##### elements

`T`[]

#### Returns

`Annoy`\<`T`\>

***

### search()

> **search**(`query`, `k?`): `object`[]

Defined in: [knn/Annoy.js:279](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L279)

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

Defined in: [knn/Annoy.js:398](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L398)

#### Parameters

##### i

`number`

##### k?

`number` = `5`

#### Returns

`object`[]

#### Overrides

`KNN.search_by_index`

***

### search\_index()

> **search\_index**(`i`, `k?`): `object`[]

Defined in: [knn/Annoy.js:410](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/Annoy.js#L410)

Alias for search_by_index for backward compatibility.

#### Parameters

##### i

`number`

Index of the query element

##### k?

`number` = `5`

Number of nearest neighbors to return

#### Returns

`object`[]
