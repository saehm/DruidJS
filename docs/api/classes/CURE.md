[@saehrimnir/druidjs](../globals.md) / CURE

# Class: CURE

Defined in: clustering/CURE.js:17

CURE (Clustering Using REpresentatives)

An efficient clustering algorithm for large databases that is robust to outliers
and identifies clusters with non-spherical shapes and wide variances in size.

## Extends

- `Clustering`

## Constructors

### Constructor

> **new CURE**(`points`, `parameters?`): `CURE`

Defined in: clustering/CURE.js:36

#### Parameters

##### points

[`InputType`](../type-aliases/InputType.md)

##### parameters?

`Partial`\<[`ParametersCURE`](../interfaces/ParametersCURE.md)\> = `{}`

#### Returns

`CURE`

#### Overrides

`Clustering.constructor`

## Properties

### \_cluster_ids

> **\_cluster_ids**: `number`[] = `[]`

Defined in: clustering/CURE.js:30

---

### \_D

> **\_D**: `number`

Defined in: clustering/Clustering.js:19

#### Inherited from

`Clustering._D`

---

### \_K

> **\_K**: `number`

Defined in: clustering/CURE.js:19

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

### \_num_representatives

> **\_num_representatives**: `number`

Defined in: clustering/CURE.js:21

---

### \_parameters

> **\_parameters**: [`ParametersCURE`](../interfaces/ParametersCURE.md)

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

### \_shrink_factor

> **\_shrink_factor**: `number`

Defined in: clustering/CURE.js:23

## Methods

### get_cluster_list()

> **get_cluster_list**(): `number`[]

Defined in: clustering/CURE.js:250

#### Returns

`number`[]

#### Overrides

`Clustering.get_cluster_list`

---

### get_clusters()

> **get_clusters**(): `number`[][]

Defined in: clustering/CURE.js:243

#### Returns

`number`[][]

#### Overrides

`Clustering.get_clusters`
