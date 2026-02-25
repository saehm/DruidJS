[@saehrimnir/druidjs](../globals.md) / HierarchicalClustering

# Class: HierarchicalClustering

Defined in: [clustering/Hierarchical_Clustering.js:18](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L18)

Hierarchical Clustering

A bottom-up approach (agglomerative) to clustering that builds a tree of clusters (dendrogram).
Supports different linkage criteria: single, complete, and average.

## Extends

- `Clustering`

## Constructors

### Constructor

> **new HierarchicalClustering**(`points`, `parameters?`): `HierarchicalClustering`

Defined in: [clustering/Hierarchical_Clustering.js:26](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L26)

#### Parameters

##### points

[`InputType`](../type-aliases/InputType.md)

Data or distance matrix if metric is 'precomputed'

##### parameters?

`Partial`\<[`ParametersHierarchicalClustering`](../interfaces/ParametersHierarchicalClustering.md)\> = `{}`

#### Returns

`HierarchicalClustering`

#### Overrides

`Clustering.constructor`

## Properties

### \_c_size

> **\_c_size**: `Uint16Array`\<`ArrayBuffer`\>

Defined in: [clustering/Hierarchical_Clustering.js:85](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L85)

---

### \_clusters

> **\_clusters**: `any`[]

Defined in: [clustering/Hierarchical_Clustering.js:83](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L83)

---

### \_D

> **\_D**: `number`

Defined in: clustering/Clustering.js:19

#### Inherited from

`Clustering._D`

---

### \_d_min

> **\_d_min**: `Float64Array`\<`ArrayBuffer`\>

Defined in: [clustering/Hierarchical_Clustering.js:41](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L41)

---

### \_distance_matrix

> **\_distance_matrix**: [`Matrix`](Matrix.md)

Defined in: [clustering/Hierarchical_Clustering.js:82](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L82)

---

### \_id

> **\_id**: `number`

Defined in: [clustering/Hierarchical_Clustering.js:33](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L33)

---

### \_matrix

> **\_matrix**: [`Matrix`](Matrix.md)

Defined in: clustering/Clustering.js:15

#### Inherited from

`Clustering._matrix`

---

### \_N

> **\_N**: `number`

Defined in: clustering/Clustering.js:17

#### Inherited from

`Clustering._N`

---

### \_parameters

> **\_parameters**: [`ParametersHierarchicalClustering`](../interfaces/ParametersHierarchicalClustering.md)

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

### root

> **root**: `Cluster` \| `null` = `null`

Defined in: [clustering/Hierarchical_Clustering.js:20](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L20)

## Methods

### get_cluster_list()

> **get_cluster_list**(`value`, `type?`): `number`[]

Defined in: [clustering/Hierarchical_Clustering.js:228](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L228)

#### Parameters

##### value

`number`

Value where to cut the tree.

##### type?

Type of value. Default is `"distance"`

`"distance"` | `"depth"`

#### Returns

`number`[]

- Array of clusters with the indices of the rows in given points.

#### Overrides

`Clustering.get_cluster_list`

---

### get_clusters()

> **get_clusters**(`value`, `type?`): `number`[][]

Defined in: [clustering/Hierarchical_Clustering.js:204](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L204)

#### Parameters

##### value

`number`

Value where to cut the tree.

##### type?

Type of value. Default is `"distance"`

`"distance"` | `"depth"`

#### Returns

`number`[][]

- Array of clusters with the indices of the rows in given points.

#### Overrides

`Clustering.get_clusters`

---

### get_clusters_raw()

> **get_clusters_raw**(`value`, `type?`): `Cluster`[][]

Defined in: [clustering/Hierarchical_Clustering.js:180](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/clustering/Hierarchical_Clustering.js#L180)

#### Parameters

##### value

`number`

Value where to cut the tree.

##### type?

Type of value. Default is `"distance"`

`"distance"` | `"depth"`

#### Returns

`Cluster`[][]

- Array of clusters with the indices of the rows in given points.
