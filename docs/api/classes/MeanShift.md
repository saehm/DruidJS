[@saehrimnir/druidjs](../globals.md) / MeanShift

# Class: MeanShift

Defined in: clustering/MeanShift.js:18

Mean Shift Clustering

A non-parametric clustering technique that does not require prior knowledge of the
number of clusters. It identifies centers of density in the data.

## Extends

- `Clustering`

## Constructors

### Constructor

> **new MeanShift**(`points`, `parameters?`): `MeanShift`

Defined in: clustering/MeanShift.js:38

#### Parameters

##### points

[`InputType`](../type-aliases/InputType.md)

##### parameters?

`Partial`\<[`ParametersMeanShift`](../interfaces/ParametersMeanShift.md)\> = `{}`

#### Returns

`MeanShift`

#### Overrides

`Clustering.constructor`

## Properties

### \_bandwidth

> **\_bandwidth**: `number`

Defined in: clustering/MeanShift.js:20

---

### \_cluster_list

> **\_cluster_list**: `number`[][] \| `undefined`

Defined in: clustering/MeanShift.js:32

---

### \_clusters

> **\_clusters**: `number`[] \| `undefined`

Defined in: clustering/MeanShift.js:30

---

### \_D

> **\_D**: `number`

Defined in: clustering/Clustering.js:19

#### Inherited from

`Clustering._D`

---

### \_kernel()

> **\_kernel**: (`dist`) => `number`

Defined in: clustering/MeanShift.js:26

#### Parameters

##### dist

`number`

#### Returns

`number`

---

### \_matrix

> **\_matrix**: [`Matrix`](Matrix.md)

Defined in: clustering/Clustering.js:15

#### Inherited from

`Clustering._matrix`

---

### \_max_iter

> **\_max_iter**: `number`

Defined in: clustering/MeanShift.js:22

---

### \_N

> **\_N**: `number`

Defined in: clustering/Clustering.js:17

#### Inherited from

`Clustering._N`

---

### \_parameters

> **\_parameters**: [`ParametersMeanShift`](../interfaces/ParametersMeanShift.md)

Defined in: clustering/Clustering.js:13

#### Inherited from

`Clustering._parameters`

---

### \_points

> **\_points**: [`Matrix`](Matrix.md)

Defined in: clustering/MeanShift.js:28

#### Overrides

`Clustering._points`

---

### \_tolerance

> **\_tolerance**: `number`

Defined in: clustering/MeanShift.js:24

## Methods

### \_assign_clusters()

> **\_assign_clusters**(): `void`

Defined in: clustering/MeanShift.js:156

#### Returns

`void`

---

### \_compute_bandwidth()

> **\_compute_bandwidth**(`matrix`): `number`

Defined in: clustering/MeanShift.js:76

#### Parameters

##### matrix

[`Matrix`](Matrix.md)

#### Returns

`number`

---

### \_kernel_weight()

> **\_kernel_weight**(`dist`): `number`

Defined in: clustering/MeanShift.js:99

#### Parameters

##### dist

`number`

#### Returns

`number`

---

### \_mean_shift()

> **\_mean_shift**(): `void`

Defined in: clustering/MeanShift.js:104

#### Returns

`void`

---

### get_cluster_list()

> **get_cluster_list**(): `number`[]

Defined in: clustering/MeanShift.js:224

#### Returns

`number`[]

#### Overrides

`Clustering.get_cluster_list`

---

### get_clusters()

> **get_clusters**(): `number`[][]

Defined in: clustering/MeanShift.js:211

#### Returns

`number`[][]

#### Overrides

`Clustering.get_clusters`
