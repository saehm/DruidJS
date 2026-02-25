[@saehrimnir/druidjs](../globals.md) / BallTree

# Class: BallTree\<T\>

Defined in: [knn/BallTree.js:27](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/BallTree.js#L27)

Ball Tree for efficient nearest neighbor search.

A Ball Tree is a metric tree that partitions points into a nested set of
hyperspheres (balls). It is particularly effective for high-dimensional
data and supports any valid metric.

## Template

## Extends

- `KNN`

## Type Parameters

### T

`T` _extends_ `number`[] \| `Float64Array`

## Constructors

### Constructor

> **new BallTree**\<`T`\>(`elements`, `parameters?`): `BallTree`\<`T`\>

Defined in: [knn/BallTree.js:36](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/BallTree.js#L36)

Generates a BallTree with given `elements`.

#### Parameters

##### elements

`T`[]

Elements which should be added to the BallTree

##### parameters?

[`ParametersBallTree`](../interfaces/ParametersBallTree.md) = `...`

Default is `{metric: euclidean}`

#### Returns

`BallTree`\<`T`\>

#### See

- [https://en.wikipedia.org/wiki/Ball_tree](https://en.wikipedia.org/wiki/Ball_tree)
- [https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js](https://github.com/invisal/noobjs/blob/master/src/tree/BallTree.js)

#### Overrides

`KNN.constructor`

## Properties

### \_elements

> **\_elements**: `T`[]

Defined in: [knn/KNN.js:14](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L14)

#### Inherited from

`KNN._elements`

---

### \_parameters

> **\_parameters**: [`ParametersBallTree`](../interfaces/ParametersBallTree.md)

Defined in: [knn/KNN.js:16](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L16)

#### Inherited from

`KNN._parameters`

---

### \_type

> **\_type**: `"array"` \| `"typed"`

Defined in: [knn/KNN.js:18](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/KNN.js#L18)

#### Inherited from

`KNN._type`

## Accessors

### \_metric

#### Get Signature

> **get** **\_metric**(): [`Metric`](../type-aliases/Metric.md)

Defined in: [knn/BallTree.js:46](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/BallTree.js#L46)

##### Returns

[`Metric`](../type-aliases/Metric.md)

## Methods

### search()

> **search**(`t`, `k?`): `object`[]

Defined in: [knn/BallTree.js:119](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/BallTree.js#L119)

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

Defined in: [knn/BallTree.js:110](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/knn/BallTree.js#L110)

#### Parameters

##### i

`number`

##### k?

`number` = `5`

#### Returns

`object`[]

#### Overrides

`KNN.search_by_index`
