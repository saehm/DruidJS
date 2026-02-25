[@saehrimnir/druidjs](../globals.md) / XMeans

# Class: XMeans

Defined in: [clustering/XMeans.js:34](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/XMeans.js#L34)

X-Means Clustering

An extension of K-Means that automatically determines the number of clusters (K)
using the Bayesian Information Criterion (BIC).

## Extends

- `Clustering`

## Constructors

### Constructor

> **new XMeans**(`points`, `parameters?`): `XMeans`

Defined in: [clustering/XMeans.js:54](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/XMeans.js#L54)

XMeans clustering algorithm that automatically determines the optimal number of clusters.

X-Means extends K-Means by starting with a minimum number of clusters and iteratively
splitting clusters to improve the Bayesian Information Criterion (BIC).

Algorithm:
1. Start with K_min clusters using KMeans
2. For each cluster, try splitting it into 2 sub-clusters
3. If BIC improves after splitting, keep the split
4. Run KMeans again with all (old + new) centroids
5. Repeat until K_max is reached or no more improvements

#### Parameters

##### points

[`InputType`](../type-aliases/InputType.md)

The data points to cluster

##### parameters?

`Partial`\<[`ParametersXMeans`](../interfaces/ParametersXMeans.md)\> = `{}`

Configuration parameters

#### Returns

`XMeans`

#### See

 - [https://www.cs.cmu.edu/~dpelleg/download/xmeans.pdf](https://www.cs.cmu.edu/~dpelleg/download/xmeans.pdf)
 - [https://github.com/annoviko/pyclustering/blob/master/pyclustering/cluster/xmeans.py](https://github.com/annoviko/pyclustering/blob/master/pyclustering/cluster/xmeans.py)
 - [https://github.com/haifengl/smile/blob/master/core/src/main/java/smile/clustering/XMeans.java](https://github.com/haifengl/smile/blob/master/core/src/main/java/smile/clustering/XMeans.java)

#### Overrides

`Clustering.constructor`

## Properties

### \_best\_kmeans

> **\_best\_kmeans**: [`KMeans`](KMeans.md) \| `null`

Defined in: [clustering/XMeans.js:67](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/XMeans.js#L67)

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

> **\_parameters**: [`ParametersXMeans`](../interfaces/ParametersXMeans.md)

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

Defined in: [clustering/XMeans.js:64](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/XMeans.js#L64)

## Accessors

### centroids

#### Get Signature

> **get** **centroids**(): `Float64Array`\<`ArrayBufferLike`\>[]

Defined in: [clustering/XMeans.js:344](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/XMeans.js#L344)

Get the final centroids

##### Returns

`Float64Array`\<`ArrayBufferLike`\>[]

Array of centroids

***

### k

#### Get Signature

> **get** **k**(): `number`

Defined in: [clustering/XMeans.js:356](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/XMeans.js#L356)

Get the optimal number of clusters found

##### Returns

`number`

The number of clusters

## Methods

### get\_cluster\_list()

> **get\_cluster\_list**(): `number`[]

Defined in: [clustering/XMeans.js:332](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/XMeans.js#L332)

#### Returns

`number`[]

The cluster list

#### Overrides

`Clustering.get_cluster_list`

***

### get\_clusters()

> **get\_clusters**(): `number`[][]

Defined in: [clustering/XMeans.js:324](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/XMeans.js#L324)

Get the computed clusters

#### Returns

`number`[][]

Array of clusters, each containing indices of points

#### Overrides

`Clustering.get_clusters`
