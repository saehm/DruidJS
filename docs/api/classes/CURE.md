[@saehrimnir/druidjs](../globals.md) / CURE

# Class: CURE

Defined in: [clustering/CURE.js:17](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/CURE.js#L17)

CURE (Clustering Using REpresentatives)

An efficient clustering algorithm for large databases that is robust to outliers
and identifies clusters with non-spherical shapes and wide variances in size.

## Extends

- `Clustering`

## Constructors

### Constructor

> **new CURE**(`points`, `parameters?`): `CURE`

Defined in: [clustering/CURE.js:36](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/CURE.js#L36)

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

### \_cluster\_ids

> **\_cluster\_ids**: `number`[] = `[]`

Defined in: [clustering/CURE.js:30](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/CURE.js#L30)

***

### \_D

> **\_D**: `number`

Defined in: [clustering/Clustering.js:19](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L19)

#### Inherited from

`Clustering._D`

***

### \_K

> **\_K**: `number`

Defined in: [clustering/CURE.js:19](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/CURE.js#L19)

***

### \_matrix

> **\_matrix**: [`Matrix`](Matrix.md)

Defined in: [clustering/Clustering.js:15](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L15)

#### Inherited from

`Clustering._matrix`

***

### \_N

> **\_N**: `number`

Defined in: [clustering/Clustering.js:17](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L17)

#### Inherited from

`Clustering._N`

***

### \_num\_representatives

> **\_num\_representatives**: `number`

Defined in: [clustering/CURE.js:21](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/CURE.js#L21)

***

### \_parameters

> **\_parameters**: [`ParametersCURE`](../interfaces/ParametersCURE.md)

Defined in: [clustering/Clustering.js:13](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L13)

#### Inherited from

`Clustering._parameters`

***

### \_points

> **\_points**: [`InputType`](../type-aliases/InputType.md)

Defined in: [clustering/Clustering.js:11](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L11)

#### Inherited from

`Clustering._points`

***

### \_shrink\_factor

> **\_shrink\_factor**: `number`

Defined in: [clustering/CURE.js:23](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/CURE.js#L23)

## Methods

### get\_cluster\_list()

> **get\_cluster\_list**(): `number`[]

Defined in: [clustering/CURE.js:250](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/CURE.js#L250)

#### Returns

`number`[]

#### Overrides

`Clustering.get_cluster_list`

***

### get\_clusters()

> **get\_clusters**(): `number`[][]

Defined in: [clustering/CURE.js:243](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/CURE.js#L243)

#### Returns

`number`[][]

#### Overrides

`Clustering.get_clusters`
