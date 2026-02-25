[@saehrimnir/druidjs](../globals.md) / TriMap

# Class: TriMap\<T\>

Defined in: [dimred/TriMap.js:24](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L24)

TriMap

A dimensionality reduction technique that preserves both local and global
structure using triplets. It is designed to be a more robust alternative
to t-SNE and UMAP.

## Template

## Extends

- `DR`

## Type Parameters

### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new TriMap**\<`T`\>(`X`, `parameters?`): `TriMap`\<`T`\>

Defined in: [dimred/TriMap.js:31](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L31)

#### Parameters

##### X

`T`

The high-dimensional data.

##### parameters?

`Partial`\<[`ParametersTriMap`](../interfaces/ParametersTriMap.md)\>

Object containing parameterization of the DR method.

#### Returns

`TriMap`\<`T`\>

#### See

 - [https://arxiv.org/pdf/1910.00204v1.pdf](https://arxiv.org/pdf/1910.00204v1.pdf)
 - [https://github.com/eamid/trimap](https://github.com/eamid/trimap)

#### Overrides

`DR.constructor`

## Properties

### \_\_input

> **\_\_input**: `T`

Defined in: [dimred/DR.js:38](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L38)

#### Inherited from

`DR.__input`

***

### \_D

> **\_D**: `number`

Defined in: [dimred/DR.js:20](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L20)

#### Inherited from

`DR._D`

***

### \_is\_initialized

> **\_is\_initialized**: `boolean`

Defined in: [dimred/DR.js:26](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L26)

#### Inherited from

`DR._is_initialized`

***

### \_N

> **\_N**: `number`

Defined in: [dimred/DR.js:22](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L22)

#### Inherited from

`DR._N`

***

### \_parameters

> **\_parameters**: [`ParametersTriMap`](../interfaces/ParametersTriMap.md)

Defined in: [dimred/DR.js:41](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L41)

#### Inherited from

`DR._parameters`

***

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: [dimred/DR.js:24](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L24)

#### Inherited from

`DR._randomizer`

***

### \_type

> **\_type**: `"array"` \| `"matrix"` \| `"typed"`

Defined in: [dimred/DR.js:46](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L46)

#### Inherited from

`DR._type`

***

### C

> **C**: `number` \| `undefined`

Defined in: [dimred/TriMap.js:68](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L68)

***

### gain

> **gain**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/TriMap.js:70](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L70)

***

### knn

> **knn**: `KNN`\<`number`[] \| `Float64Array`\<`ArrayBufferLike`\>, `any`\> \| `undefined`

Defined in: [dimred/TriMap.js:63](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L63)

***

### lr

> **lr**: `number` \| `undefined`

Defined in: [dimred/TriMap.js:67](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L67)

***

### n\_inliers

> **n\_inliers**: `number` \| `undefined`

Defined in: [dimred/TriMap.js:59](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L59)

***

### n\_outliers

> **n\_outliers**: `number` \| `undefined`

Defined in: [dimred/TriMap.js:60](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L60)

***

### n\_random

> **n\_random**: `number` \| `undefined`

Defined in: [dimred/TriMap.js:61](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L61)

***

### triplets

> **triplets**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/TriMap.js:65](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L65)

***

### vel

> **vel**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/TriMap.js:69](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L69)

***

### weights

> **weights**: `Float64Array`\<`ArrayBuffer`\> \| `undefined`

Defined in: [dimred/TriMap.js:66](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L66)

***

### X

> **X**: [`Matrix`](Matrix.md)

Defined in: [dimred/DR.js:48](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L48)

#### Inherited from

`DR.X`

***

### Y

> **Y**: [`Matrix`](Matrix.md)

Defined in: [dimred/DR.js:50](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L50)

#### Inherited from

`DR.Y`

## Accessors

### projection

#### Get Signature

> **get** **projection**(): `T`

Defined in: [dimred/DR.js:211](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L211)

##### Returns

`T`

The projection in the type of input `X`.

#### Inherited from

`DR.projection`

## Methods

### \_generate\_triplets()

> **\_generate\_triplets**(`n_inliers`, `n_outliers`, `n_random`): `object`

Defined in: [dimred/TriMap.js:81](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L81)

Generates [n\_inliers](#generate-triplets) x [n\_outliers](#generate-triplets) x [n\_random](#generate-triplets) triplets.

#### Parameters

##### n\_inliers

`number`

##### n\_outliers

`number`

##### n\_random

`number`

#### Returns

`object`

##### triplets

> **triplets**: [`Matrix`](Matrix.md)

##### weights

> **weights**: `Float64Array`\<`ArrayBuffer`\>

***

### \_grad()

> **\_grad**(`Y`): `object`

Defined in: [dimred/TriMap.js:299](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L299)

Computes the gradient for updating the embedding.

#### Parameters

##### Y

[`Matrix`](Matrix.md)

The embedding

#### Returns

`object`

##### grad

> **grad**: [`Matrix`](Matrix.md)

##### loss

> **loss**: `number`

##### n\_viol

> **n\_viol**: `number`

***

### check\_init()

> **check\_init**(): `DR`\<`T`, [`ParametersTriMap`](../interfaces/ParametersTriMap.md)\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersTriMap`](../interfaces/ParametersTriMap.md)\>

#### Inherited from

`DR.check_init`

***

### generator()

> **generator**(`max_iteration?`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/TriMap.js:373](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L373)

#### Parameters

##### max\_iteration?

`number` = `800`

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

***

### init()

> **init**(`pca?`, `knn?`): `TriMap`\<`T`\>

Defined in: [dimred/TriMap.js:52](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L52)

#### Parameters

##### pca?

Initial Embedding (if null then PCA gets used). Default is `null`

[`Matrix`](Matrix.md) | `null`

##### knn?

KNN Object (if null then BallTree gets used). Default is `null`

`KNN`\<`number`[] \| `Float64Array`\<`ArrayBufferLike`\>, `any`\> | `null`

#### Returns

`TriMap`\<`T`\>

#### Overrides

`DR.init`

***

### parameter()

#### Call Signature

> **parameter**(): [`ParametersTriMap`](../interfaces/ParametersTriMap.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersTriMap`](../interfaces/ParametersTriMap.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersTriMap`](../interfaces/ParametersTriMap.md)\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersTriMap`](../interfaces/ParametersTriMap.md)

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersTriMap`](../interfaces/ParametersTriMap.md)\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `TriMap`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersTriMap`](../interfaces/ParametersTriMap.md)

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersTriMap`](../interfaces/ParametersTriMap.md)\[`K`\]

Value of the parameter to set.

##### Returns

`TriMap`\<`T`\>

##### Inherited from

`DR.parameter`

***

### transform()

> **transform**(`max_iteration?`): `T`

Defined in: [dimred/TriMap.js:361](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L361)

#### Parameters

##### max\_iteration?

`number` = `800`

#### Returns

`T`

#### Overrides

`DR.transform`

***

### transform\_async()

> **transform\_async**(...`args`): `Promise`\<`T`\>

Defined in: [dimred/DR.js:233](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L233)

Computes the projection.

#### Parameters

##### args

...`unknown`[]

Arguments the transform method of the respective DR method takes.

#### Returns

`Promise`\<`T`\>

The dimensionality reduced dataset.

#### Inherited from

`DR.transform_async`

***

### generator()

> `static` **generator**\<`T`\>(`X`, `parameters?`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/TriMap.js:449](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L449)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersTriMap`](../interfaces/ParametersTriMap.md)\>

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

***

### transform()

> `static` **transform**\<`T`\>(`X`, `parameters?`): `T`

Defined in: [dimred/TriMap.js:438](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L438)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersTriMap`](../interfaces/ParametersTriMap.md)\>

#### Returns

`T`

#### Overrides

`DR.transform`

***

### transform\_async()

> `static` **transform\_async**\<`T`\>(`X`, `parameters?`): `Promise`\<`T`\>

Defined in: [dimred/TriMap.js:461](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TriMap.js#L461)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersTriMap`](../interfaces/ParametersTriMap.md)\>

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
