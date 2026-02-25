[@saehrimnir/druidjs](../globals.md) / KMeans

# Class: KMeans

Defined in: [clustering/KMeans.js:29](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L29)

K-Means Clustering

A popular clustering algorithm that partitions data into K clusters where each point
belongs to the cluster with the nearest mean (centroid).

## See

[KMedoids](KMedoids.md) for a more robust alternative

## Example

```ts
import * as druid from "@saehrimnir/druidjs";

const points = [
  [1, 1],
  [1.5, 1.5],
  [5, 5],
  [5.5, 5.5],
];
const kmeans = new druid.KMeans(points, { K: 2 });

const clusters = kmeans.get_cluster_list(); // [0, 0, 1, 1]
const centroids = kmeans.centroids; // center points
```

## Extends

- `Clustering`

## Constructors

### Constructor

> **new KMeans**(`points`, `parameters?`): `KMeans`

Defined in: [clustering/KMeans.js:34](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L34)

#### Parameters

##### points

[`InputType`](../type-aliases/InputType.md)

##### parameters?

`Partial`\<[`ParametersKMeans`](../interfaces/ParametersKMeans.md)\> = `{}`

#### Returns

`KMeans`

#### Overrides

`Clustering.constructor`

## Properties

### \_cluster_centroids

> **\_cluster_centroids**: `Float64Array`\<`ArrayBufferLike`\>[]

Defined in: [clustering/KMeans.js:60](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L60)

---

### \_clusters

> **\_clusters**: `number`[]

Defined in: [clustering/KMeans.js:58](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L58)

---

### \_D

> **\_D**: `number`

Defined in: [clustering/KMeans.js:52](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L52)

#### Inherited from

`Clustering._D`

---

### \_K

> **\_K**: `number`

Defined in: [clustering/KMeans.js:54](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L54)

---

### \_matrix

> **\_matrix**: [`Matrix`](Matrix.md)

Defined in: [clustering/KMeans.js:45](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L45)

#### Inherited from

`Clustering._matrix`

---

### \_N

> **\_N**: `number`

Defined in: [clustering/KMeans.js:51](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L51)

#### Inherited from

`Clustering._N`

---

### \_parameters

> **\_parameters**: [`ParametersKMeans`](../interfaces/ParametersKMeans.md)

Defined in: clustering/Clustering.js:13

#### Inherited from

`Clustering._parameters`

---

### \_points

> **\_points**: [`InputType`](../type-aliases/InputType.md)

Defined in: clustering/Clustering.js:11

#### Inherited from

`Clustering._points`

---

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: [clustering/KMeans.js:55](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L55)

## Accessors

### centroids

#### Get Signature

> **get** **centroids**(): `Float64Array`\<`ArrayBufferLike`\>[]

Defined in: [clustering/KMeans.js:84](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L84)

##### Returns

`Float64Array`\<`ArrayBufferLike`\>[]

The cluster centroids

---

### k

#### Get Signature

> **get** **k**(): `number`

Defined in: [clustering/KMeans.js:79](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L79)

##### Returns

`number`

The number of clusters

## Methods

### get_cluster_list()

> **get_cluster_list**(): `number`[]

Defined in: [clustering/KMeans.js:89](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L89)

#### Returns

`number`[]

The cluster list

#### Overrides

`Clustering.get_cluster_list`

---

### get_clusters()

> **get_clusters**(): `number`[][]

Defined in: [clustering/KMeans.js:94](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/KMeans.js#L94)

#### Returns

`number`[][]

An Array of clusters with the indices of the points.

#### Overrides

`Clustering.get_clusters`
