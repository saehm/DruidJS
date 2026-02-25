[@saehrimnir/druidjs](../globals.md) / KMedoids

# Class: KMedoids

Defined in: [clustering/KMedoids.js:20](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L20)

K-Medoids (PAM - Partitioning Around Medoids)

A robust clustering algorithm similar to K-Means, but uses actual data points (medoids)
as cluster centers and can work with any distance metric.

## See

[KMeans](KMeans.md) for a faster but less robust alternative

## Extends

- `Clustering`

## Constructors

### Constructor

> **new KMedoids**(`points`, `parameters?`): `KMedoids`

Defined in: [clustering/KMedoids.js:26](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L26)

#### Parameters

##### points

[`InputType`](../type-aliases/InputType.md)

Data matrix

##### parameters?

`Partial`\<[`ParametersKMedoids`](../interfaces/ParametersKMedoids.md)\> = `{}`

#### Returns

`KMedoids`

#### See

[https://link.springer.com/chapter/10.1007/978-3-030-32047-8\_16](https://link.springer.com/chapter/10.1007/978-3-030-32047-8_16) Faster k-Medoids Clustering: Improving the PAM, CLARA, and CLARANS Algorithms

#### Overrides

`Clustering.constructor`

## Properties

### \_A

> **\_A**: `Float64Array`\<`ArrayBufferLike`\>[]

Defined in: [clustering/KMedoids.js:28](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L28)

***

### \_cluster\_medoids

> **\_cluster\_medoids**: `number`[]

Defined in: [clustering/KMedoids.js:39](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L39)

***

### \_clusters

> **\_clusters**: `any`[]

Defined in: [clustering/KMedoids.js:38](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L38)

***

### \_D

> **\_D**: `number`

Defined in: [clustering/Clustering.js:19](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L19)

#### Inherited from

`Clustering._D`

***

### \_distance\_matrix

> **\_distance\_matrix**: [`Matrix`](Matrix.md)

Defined in: [clustering/KMedoids.js:32](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L32)

***

### \_is\_initialized

> **\_is\_initialized**: `boolean`

Defined in: [clustering/KMedoids.js:40](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L40)

***

### \_matrix

> **\_matrix**: [`Matrix`](Matrix.md)

Defined in: [clustering/Clustering.js:15](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L15)

#### Inherited from

`Clustering._matrix`

***

### \_max\_iter

> **\_max\_iter**: `number`

Defined in: [clustering/KMedoids.js:31](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L31)

***

### \_N

> **\_N**: `number`

Defined in: [clustering/Clustering.js:17](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/Clustering.js#L17)

#### Inherited from

`Clustering._N`

***

### \_parameters

> **\_parameters**: [`ParametersKMedoids`](../interfaces/ParametersKMedoids.md)

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

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: [clustering/KMedoids.js:37](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L37)

## Accessors

### k

#### Get Signature

> **get** **k**(): `number`

Defined in: [clustering/KMedoids.js:71](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L71)

##### Returns

`number`

***

### medoids

#### Get Signature

> **get** **medoids**(): `number`[]

Defined in: [clustering/KMedoids.js:76](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L76)

##### Returns

`number`[]

## Methods

### \_get\_distance()

> **\_get\_distance**(`i`, `j`, `x_i?`, `x_j?`): `number`

Defined in: [clustering/KMedoids.js:231](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L231)

#### Parameters

##### i

`number`

##### j

`number`

##### x\_i?

`Float64Array`\<`ArrayBufferLike`\> | `null`

##### x\_j?

`Float64Array`\<`ArrayBufferLike`\> | `null`

#### Returns

`number`

***

### \_get\_random\_medoids()

> **\_get\_random\_medoids**(`K`): `number`[]

Defined in: [clustering/KMedoids.js:322](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L322)

Algorithm 3. FastPAM LAB: Linear Approximate BUILD initialization.

#### Parameters

##### K

`number`

Number of clusters

#### Returns

`number`[]

***

### \_iteration()

> **\_iteration**(): `boolean`

Defined in: [clustering/KMedoids.js:161](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L161)

FastPAM1: One best swap per iteration

#### Returns

`boolean`

***

### \_nearest\_medoid()

> **\_nearest\_medoid**(`x_j`, `j`): `object`

Defined in: [clustering/KMedoids.js:251](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L251)

#### Parameters

##### x\_j

`Float64Array`\<`ArrayBufferLike`\>

##### j

`number`

#### Returns

`object`

##### distance\_nearest

> **distance\_nearest**: `number` = `d_n`

##### distance\_second

> **distance\_second**: `number` = `d_s`

##### index\_nearest

> **index\_nearest**: `number` = `n`

##### index\_second

> **index\_second**: `number` = `s`

***

### \_update\_clusters()

> **\_update\_clusters**(): `void`

Defined in: [clustering/KMedoids.js:287](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L287)

#### Returns

`void`

***

### generator()

> **generator**(): `AsyncGenerator`\<`number`[][], `void`, `unknown`\>

Defined in: [clustering/KMedoids.js:89](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L89)

#### Returns

`AsyncGenerator`\<`number`[][], `void`, `unknown`\>

***

### get\_cluster\_list()

> **get\_cluster\_list**(): `number`[]

Defined in: [clustering/KMedoids.js:44](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L44)

#### Returns

`number`[]

The cluster list

#### Overrides

`Clustering.get_cluster_list`

***

### get\_clusters()

> **get\_clusters**(): `number`[][]

Defined in: [clustering/KMedoids.js:52](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L52)

#### Returns

`number`[][]

- Array of clusters with the indices of the rows in given points.

#### Overrides

`Clustering.get_clusters`

***

### get\_medoids()

> **get\_medoids**(): `number`[]

Defined in: [clustering/KMedoids.js:81](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L81)

#### Returns

`number`[]

***

### init()

> **init**(`K`, `cluster_medoids`): `KMedoids`

Defined in: [clustering/KMedoids.js:301](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/KMedoids.js#L301)

Computes `K` clusters out of the `matrix`.

#### Parameters

##### K

`number`

Number of clusters.

##### cluster\_medoids

`number`[]

#### Returns

`KMedoids`
