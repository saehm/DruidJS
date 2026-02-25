[@saehrimnir/druidjs](../globals.md) / NNDescent

# Class: NNDescent\<T\>

Defined in: [knn/NNDescent.js:38](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/NNDescent.js#L38)

NN-Descent

An efficient graph-based approximate nearest neighbor search algorithm.
It works by iteratively improving a neighbor graph using the fact that
"neighbors of neighbors are likely to be neighbors".

## Template

## See

[Paper](http://www.cs.princeton.edu/cass/papers/www11.pdf|NN-Descent)

## Extends

- `KNN`

## Type Parameters

### T

`T` *extends* `number`[] \| `Float64Array`

## Constructors

### Constructor

> **new NNDescent**\<`T`\>(`elements`, `parameters?`): `NNDescent`\<`T`\>

Defined in: [knn/NNDescent.js:55](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/NNDescent.js#L55)

#### Parameters

##### elements

`T`[]

Called V in paper.

##### parameters?

`Partial`\<[`ParametersNNDescent`](../interfaces/ParametersNNDescent.md)\> = `{}`

#### Returns

`NNDescent`\<`T`\>

#### See

[http://www.cs.princeton.edu/cass/papers/www11.pdf](http://www.cs.princeton.edu/cass/papers/www11.pdf)

#### Overrides

`KNN.constructor`

## Properties

### \_elements

> **\_elements**: `T`[]

Defined in: [knn/KNN.js:14](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L14)

#### Inherited from

`KNN._elements`

***

### \_N

> **\_N**: `number`

Defined in: [knn/NNDescent.js:62](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/NNDescent.js#L62)

***

### \_nndescent\_elements

> **\_nndescent\_elements**: `object`[]

Defined in: [knn/NNDescent.js:66](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/NNDescent.js#L66)

#### flag

> **flag**: `boolean` = `true`

#### index

> **index**: `number` = `i`

#### value

> **value**: `T` = `e`

***

### \_parameters

> **\_parameters**: [`ParametersNNDescent`](../interfaces/ParametersNNDescent.md)

Defined in: [knn/KNN.js:16](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L16)

#### Inherited from

`KNN._parameters`

***

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: [knn/NNDescent.js:63](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/NNDescent.js#L63)

***

### \_sample\_size

> **\_sample\_size**: `number`

Defined in: [knn/NNDescent.js:64](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/NNDescent.js#L64)

***

### \_type

> **\_type**: `"array"` \| `"typed"`

Defined in: [knn/KNN.js:18](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/KNN.js#L18)

#### Inherited from

`KNN._type`

## Methods

### add()

> **add**(`elements`): `NNDescent`\<`T`\>

Defined in: [knn/NNDescent.js:152](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/NNDescent.js#L152)

#### Parameters

##### elements

`T`[]

#### Returns

`NNDescent`\<`T`\>

***

### search()

> **search**(`x`, `k?`): `object`[]

Defined in: [knn/NNDescent.js:236](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/NNDescent.js#L236)

#### Parameters

##### x

`T`

##### k?

`number` = `5`

Default is `5`

#### Returns

`object`[]

#### Overrides

`KNN.search`

***

### search\_by\_index()

> **search\_by\_index**(`i`, `k?`): `object`[]

Defined in: [knn/NNDescent.js:325](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/knn/NNDescent.js#L325)

#### Parameters

##### i

`number`

##### k?

`number` = `5`

Default is `5`

#### Returns

`object`[]

#### Overrides

`KNN.search_by_index`
