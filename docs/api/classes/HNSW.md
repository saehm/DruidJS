[@saehrimnir/druidjs](../globals.md) / HNSW

# Class: HNSW\<T\>

Defined in: [knn/HNSW.js:61](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L61)

Hierarchical Navigable Small World (HNSW) graph for approximate nearest neighbor search.

HNSW builds a multi-layer graph structure where each layer is a navigable small world graph.
The top layers serve as "highways" for fast traversal, while lower layers provide accuracy.
Each element is assigned to a random level, allowing logarithmic search complexity.

Key parameters:
- `m`: Controls the number of connections per element (affects accuracy/memory)
- `ef_construction`: Controls the quality of the graph during construction (higher = better but slower)
- `ef`: Controls the quality of search (higher = better recall but slower)

Based on:
- "Efficient and robust approximate nearest neighbor search using Hierarchical Navigable Small World graphs"
  by Malkov & Yashunin (2016)
- "Approximate Nearest Neighbor Search on High Dimensional Data"
  by Li et al. (2019)

## Template

## Example

```ts
import * as druid from "@saehrimnir/druidjs";

const points = [[1, 2], [3, 4], [5, 6], [7, 8]];
const hnsw = new druid.HNSW(points, {
    metric: druid.euclidean,
    m: 16,
    ef_construction: 200
});

const query = [2, 3];
const neighbors = hnsw.search(query, 2);
// [{ element: [1, 2], index: 0, distance: 1.41 }, ...]
```

## Extends

- `KNN`

## Type Parameters

### T

`T` *extends* `number`[] \| `Float64Array`

## Constructors

### Constructor

> **new HNSW**\<`T`\>(`points`, `parameters?`): `HNSW`\<`T`\>

Defined in: [knn/HNSW.js:68](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L68)

Creates a new HNSW index.

#### Parameters

##### points

`T`[]

Initial points to add to the index

##### parameters?

[`ParametersHNSW`](../interfaces/ParametersHNSW.md) = `...`

Configuration parameters

#### Returns

`HNSW`\<`T`\>

#### Overrides

`KNN.constructor`

## Properties

### \_ef

> **\_ef**: `number`

Defined in: [knn/HNSW.js:142](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L142)

***

### \_ef\_construction

> **\_ef\_construction**: `number`

Defined in: [knn/HNSW.js:135](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L135)

***

### \_elements

> **\_elements**: `T`[]

Defined in: [knn/KNN.js:14](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L14)

#### Inherited from

`KNN._elements`

***

### \_ep

> **\_ep**: `number`[] \| `null`

Defined in: [knn/HNSW.js:161](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L161)

***

### \_L

> **\_L**: `number`

Defined in: [knn/HNSW.js:158](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L158)

***

### \_m

> **\_m**: `number`

Defined in: [knn/HNSW.js:128](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L128)

***

### \_m0

> **\_m0**: `number`

Defined in: [knn/HNSW.js:149](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L149)

***

### \_metric

> **\_metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [knn/HNSW.js:108](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L108)

***

### \_mL

> **\_mL**: `number`

Defined in: [knn/HNSW.js:152](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L152)

***

### \_next\_index

> **\_next\_index**: `number`

Defined in: [knn/HNSW.js:120](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L120)

***

### \_parameters

> **\_parameters**: [`ParametersHNSW`](../interfaces/ParametersHNSW.md)

Defined in: [knn/KNN.js:16](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L16)

#### Inherited from

`KNN._parameters`

***

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: [knn/HNSW.js:155](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L155)

***

### \_select

> **\_select**: `Function`

Defined in: [knn/HNSW.js:111](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L111)

***

### \_type

> **\_type**: `"array"` \| `"typed"`

Defined in: [knn/KNN.js:18](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L18)

#### Inherited from

`KNN._type`

## Accessors

### num\_layers

#### Get Signature

> **get** **num\_layers**(): `number`

Defined in: [knn/HNSW.js:710](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L710)

Get the number of layers in the graph.

##### Returns

`number`

Number of layers

***

### size

#### Get Signature

> **get** **size**(): `number`

Defined in: [knn/HNSW.js:701](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L701)

Get the number of elements in the index.

##### Returns

`number`

Number of elements

## Methods

### add()

> **add**(`new_elements`): `HNSW`\<`T`\>

Defined in: [knn/HNSW.js:185](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L185)

Add multiple elements to the index.

#### Parameters

##### new\_elements

`T`[]

Elements to add

#### Returns

`HNSW`\<`T`\>

This instance for chaining

***

### addOne()

> **addOne**(`element`): `HNSW`\<`T`\>

Defined in: [knn/HNSW.js:175](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L175)

Add a single element to the index.

#### Parameters

##### element

`T`

Element to add

#### Returns

`HNSW`\<`T`\>

This instance for chaining

***

### get\_element()

> **get\_element**(`index`): `T`

Defined in: [knn/HNSW.js:720](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L720)

Get an element by its index.

#### Parameters

##### index

`number`

Element index

#### Returns

`T`

The element at the given index

***

### search()

> **search**(`q`, `K`): `Candidate`\<`T`\>[]

Defined in: [knn/HNSW.js:578](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L578)

Searches for the K nearest neighbors to a query element in the HNSW graph.

Performs a multi-layer search starting from the entry point and traversing
each layer as entry points for the next.

#### Parameters

##### q

`T`

Query element

##### K

`number`

Number of nearest neighbors to return

#### Returns

`Candidate`\<`T`\>[]

K nearest neighbors with their distances

#### Overrides

`KNN.search`

***

### search\_by\_index()

> **search\_by\_index**(`i`, `K?`): `Candidate`\<`T`\>[]

Defined in: [knn/HNSW.js:731](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L731)

Search for nearest neighbors using an element index as the query.

#### Parameters

##### i

`number`

Index of the query element

##### K?

`number` = `5`

Number of nearest neighbors to return

#### Returns

`Candidate`\<`T`\>[]

K nearest neighbors

#### Overrides

`KNN.search_by_index`

***

### search\_iter()

> **search\_iter**(`q`, `K`, `ef?`): `Generator`\<\{ `candidates`: `object`[]; `layer`: `number`; \}, `void`, `unknown`\>

Defined in: [knn/HNSW.js:660](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/HNSW.js#L660)

Iterator for searching the HNSW graph layer by layer.

Yields intermediate results at each layer for debugging or visualization.

#### Parameters

##### q

`T`

Query element

##### K

`number`

Number of nearest neighbors to return

##### ef?

Size of dynamic candidate list

`number` | `null`

#### Returns

`Generator`\<\{ `candidates`: `object`[]; `layer`: `number`; \}, `void`, `unknown`\>

#### Yields
