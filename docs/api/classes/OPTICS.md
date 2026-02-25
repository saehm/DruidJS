[@saehrimnir/druidjs](../globals.md) / OPTICS

# Class: OPTICS

Defined in: [clustering/OPTICS.js:26](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/OPTICS.js#L26)

OPTICS (Ordering Points To Identify the Clustering Structure)

A density-based clustering algorithm that extends DBSCAN. It handles clusters of varying
densities and produces a reachability plot that can be used to extract clusters.

## Extends

- `Clustering`

## Constructors

### Constructor

> **new OPTICS**(`points`, `parameters?`): `OPTICS`

Defined in: [clustering/OPTICS.js:35](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/OPTICS.js#L35)

**O**rdering **P**oints **T**o **I**dentify the **C**lustering **S**tructure.

#### Parameters

##### points

[`InputType`](../type-aliases/InputType.md)

The data.

##### parameters?

`Partial`\<[`ParametersOptics`](../interfaces/ParametersOptics.md)\> = `{}`

#### Returns

`OPTICS`

#### See

 - [https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf](https://www.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf)
 - [https://en.wikipedia.org/wiki/OPTICS\_algorithm](https://en.wikipedia.org/wiki/OPTICS_algorithm)

#### Overrides

`Clustering.constructor`

## Properties

### \_cluster\_index

> **\_cluster\_index**: `number`

Defined in: [clustering/OPTICS.js:69](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/OPTICS.js#L69)

***

### \_clusters

> **\_clusters**: `number`[][]

Defined in: [clustering/OPTICS.js:50](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/OPTICS.js#L50)

***

### \_D

> **\_D**: `number`

Defined in: [clustering/Clustering.js:19](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L19)

#### Inherited from

`Clustering._D`

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

### \_parameters

> **\_parameters**: [`ParametersOptics`](../interfaces/ParametersOptics.md)

Defined in: [clustering/Clustering.js:13](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L13)

#### Inherited from

`Clustering._parameters`

***

### \_points

> **\_points**: [`InputType`](../type-aliases/InputType.md)

Defined in: [clustering/Clustering.js:11](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L11)

#### Inherited from

`Clustering._points`

## Methods

### get\_cluster\_list()

> **get\_cluster\_list**(): `number`[]

Defined in: [clustering/OPTICS.js:206](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/OPTICS.js#L206)

#### Returns

`number`[]

Returns an array, where the ith entry defines the cluster affirmation of the ith point of
  given data. (-1 stands for outlier)

#### Overrides

`Clustering.get_cluster_list`

***

### get\_clusters()

> **get\_clusters**(): `number`[][]

Defined in: [clustering/OPTICS.js:187](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/OPTICS.js#L187)

Returns an array of clusters.

#### Returns

`number`[][]

Array of clusters with the indices of the rows in given `matrix`.

#### Overrides

`Clustering.get_clusters`
