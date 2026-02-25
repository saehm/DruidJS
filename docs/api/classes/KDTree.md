[@saehrimnir/druidjs](../globals.md) / KDTree

# Class: KDTree\<T\>

Defined in: [knn/KDTree.js:36](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KDTree.js#L36)

KD-Tree (K-dimensional Tree) for efficient nearest neighbor search.

KD-Trees partition k-dimensional space by recursively splitting along coordinate axes.
At each level, the tree splits points based on the median of the coordinate with the largest spread.
This creates a balanced binary tree structure that enables efficient O(log n) search on average.

Best suited for:
- Low to moderate dimensional data (d < 20-30)
- When exact nearest neighbors are needed
- When dimensionality is not too high

Performance degrades in high dimensions (curse of dimensionality) where approximate
methods like HNSW or LSH become more effective.

## Template

## See

[https://en.wikipedia.org/wiki/K-d\_tree](https://en.wikipedia.org/wiki/K-d_tree)

## Extends

- `KNN`

## Type Parameters

### T

`T` *extends* `number`[] \| `Float64Array`

## Constructors

### Constructor

> **new KDTree**\<`T`\>(`elements`, `parameters?`): `KDTree`\<`T`\>

Defined in: [knn/KDTree.js:43](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KDTree.js#L43)

Generates a KD-Tree with given `elements`.

#### Parameters

##### elements

`T`[]

Elements which should be added to the KD-Tree

##### parameters?

[`ParametersKDTree`](../interfaces/ParametersKDTree.md) = `...`

Default is `{metric: euclidean}`

#### Returns

`KDTree`\<`T`\>

#### Overrides

`KNN.constructor`

## Properties

### \_elements

> **\_elements**: `T`[]

Defined in: [knn/KNN.js:14](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L14)

#### Inherited from

`KNN._elements`

***

### \_parameters

> **\_parameters**: [`ParametersKDTree`](../interfaces/ParametersKDTree.md)

Defined in: [knn/KNN.js:16](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L16)

#### Inherited from

`KNN._parameters`

***

### \_type

> **\_type**: `"array"` \| `"typed"`

Defined in: [knn/KNN.js:18](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L18)

#### Inherited from

`KNN._type`

## Accessors

### \_metric

#### Get Signature

> **get** **\_metric**(): [`Metric`](../type-aliases/Metric.md)

Defined in: [knn/KDTree.js:56](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KDTree.js#L56)

##### Returns

[`Metric`](../type-aliases/Metric.md)

## Methods

### search()

> **search**(`t`, `k?`): `object`[]

Defined in: [knn/KDTree.js:106](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KDTree.js#L106)

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

***

### search\_by\_index()

> **search\_by\_index**(`i`, `k?`): `object`[]

Defined in: [knn/KDTree.js:97](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KDTree.js#L97)

#### Parameters

##### i

`number`

##### k?

`number` = `5`

#### Returns

`object`[]

#### Overrides

`KNN.search_by_index`
