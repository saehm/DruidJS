[@saehrimnir/druidjs](../globals.md) / MeanShift

# Class: MeanShift

Defined in: [clustering/MeanShift.js:18](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L18)

Mean Shift Clustering

A non-parametric clustering technique that does not require prior knowledge of the
number of clusters. It identifies centers of density in the data.

## Extends

- `Clustering`

## Constructors

### Constructor

> **new MeanShift**(`points`, `parameters?`): `MeanShift`

Defined in: [clustering/MeanShift.js:38](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L38)

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

Defined in: [clustering/MeanShift.js:20](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L20)

***

### \_cluster\_list

> **\_cluster\_list**: `number`[][] \| `undefined`

Defined in: [clustering/MeanShift.js:32](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L32)

***

### \_clusters

> **\_clusters**: `number`[] \| `undefined`

Defined in: [clustering/MeanShift.js:30](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L30)

***

### \_D

> **\_D**: `number`

Defined in: [clustering/Clustering.js:19](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L19)

#### Inherited from

`Clustering._D`

***

### \_kernel()

> **\_kernel**: (`dist`) => `number`

Defined in: [clustering/MeanShift.js:26](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L26)

#### Parameters

##### dist

`number`

#### Returns

`number`

***

### \_matrix

> **\_matrix**: [`Matrix`](Matrix.md)

Defined in: [clustering/Clustering.js:15](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L15)

#### Inherited from

`Clustering._matrix`

***

### \_max\_iter

> **\_max\_iter**: `number`

Defined in: [clustering/MeanShift.js:22](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L22)

***

### \_N

> **\_N**: `number`

Defined in: [clustering/Clustering.js:17](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L17)

#### Inherited from

`Clustering._N`

***

### \_parameters

> **\_parameters**: [`ParametersMeanShift`](../interfaces/ParametersMeanShift.md)

Defined in: [clustering/Clustering.js:13](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L13)

#### Inherited from

`Clustering._parameters`

***

### \_points

> **\_points**: [`Matrix`](Matrix.md)

Defined in: [clustering/MeanShift.js:28](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L28)

#### Overrides

`Clustering._points`

***

### \_tolerance

> **\_tolerance**: `number`

Defined in: [clustering/MeanShift.js:24](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L24)

## Methods

### \_assign\_clusters()

> **\_assign\_clusters**(): `void`

Defined in: [clustering/MeanShift.js:156](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L156)

#### Returns

`void`

***

### \_compute\_bandwidth()

> **\_compute\_bandwidth**(`matrix`): `number`

Defined in: [clustering/MeanShift.js:76](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L76)

#### Parameters

##### matrix

[`Matrix`](Matrix.md)

#### Returns

`number`

***

### \_kernel\_weight()

> **\_kernel\_weight**(`dist`): `number`

Defined in: [clustering/MeanShift.js:99](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L99)

#### Parameters

##### dist

`number`

#### Returns

`number`

***

### \_mean\_shift()

> **\_mean\_shift**(): `void`

Defined in: [clustering/MeanShift.js:104](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L104)

#### Returns

`void`

***

### get\_cluster\_list()

> **get\_cluster\_list**(): `number`[]

Defined in: [clustering/MeanShift.js:224](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L224)

#### Returns

`number`[]

#### Overrides

`Clustering.get_cluster_list`

***

### get\_clusters()

> **get\_clusters**(): `number`[][]

Defined in: [clustering/MeanShift.js:211](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/MeanShift.js#L211)

#### Returns

`number`[][]

#### Overrides

`Clustering.get_clusters`
