[@saehrimnir/druidjs](../globals.md) / Annoy

# Class: Annoy\<T\>

Defined in: knn/Annoy.js:45

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

`T` _extends_ `number`[] \| `Float64Array`

## Constructors

### Constructor

> **new Annoy**\<`T`\>(`elements`, `parameters?`): `Annoy`\<`T`\>

Defined in: knn/Annoy.js:52

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

Defined in: [knn/KNN.js:14](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L14)

#### Inherited from

`KNN._elements`

---

### \_maxPointsPerLeaf

> **\_maxPointsPerLeaf**: `number`

Defined in: knn/Annoy.js:69

---

### \_metric

> **\_metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: knn/Annoy.js:67

---

### \_numTrees

> **\_numTrees**: `number`

Defined in: knn/Annoy.js:68

---

### \_parameters

> **\_parameters**: [`ParametersAnnoy`](../interfaces/ParametersAnnoy.md)

Defined in: [knn/KNN.js:16](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L16)

#### Inherited from

`KNN._parameters`

---

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: knn/Annoy.js:71

---

### \_seed

> **\_seed**: `number`

Defined in: knn/Annoy.js:70

---

### \_type

> **\_type**: `"array"` \| `"typed"`

Defined in: [knn/KNN.js:18](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L18)

#### Inherited from

`KNN._type`

## Accessors

### num_nodes

#### Get Signature

> **get** **num_nodes**(): `number`

Defined in: knn/Annoy.js:101

Get the total number of nodes in all trees.

##### Returns

`number`

---

### num_trees

#### Get Signature

> **get** **num_trees**(): `number`

Defined in: knn/Annoy.js:93

Get the number of trees in the index.

##### Returns

`number`

## Methods

### add()

> **add**(`elements`): `Annoy`\<`T`\>

Defined in: knn/Annoy.js:124

Add elements to the Annoy index.

#### Parameters

##### elements

`T`[]

#### Returns

`Annoy`\<`T`\>

---

### search()

> **search**(`query`, `k?`): `object`[]

Defined in: knn/Annoy.js:279

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

Defined in: knn/Annoy.js:398

#### Parameters

##### i

`number`

##### k?

`number` = `5`

#### Returns

`object`[]

#### Overrides

`KNN.search_by_index`

---

### search_index()

> **search_index**(`i`, `k?`): `object`[]

Defined in: knn/Annoy.js:410

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
