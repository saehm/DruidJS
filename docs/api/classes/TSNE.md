[@saehrimnir/druidjs](../globals.md) / TSNE

# Class: TSNE\<T\>

Defined in: [dimred/TSNE.js:36](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L36)

t-SNE (t-Distributed Stochastic Neighbor Embedding)

A nonlinear dimensionality reduction technique particularly well-suited
for visualizing high-dimensional data in 2D or 3D. Preserves local
structure while revealing global patterns.

## Template

## See

 - [Paper](https://lvdmaaten.github.io/tsne/|t-SNE)
 - [UMAP](UMAP.md) for faster alternative with similar results

## Example

```ts
import * as druid from "@saehrimnir/druidjs";

const X = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]];
const tsne = new druid.TSNE(X, {
    perplexity: 30,
    epsilon: 10,
    d: 2,
    seed: 42
});

const Y = tsne.transform(500); // 500 iterations
// [[x1, y1], [x2, y2], [x3, y3]]
```

## Extends

- `DR`

## Type Parameters

### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new TSNE**\<`T`\>(`X`, `parameters?`): `TSNE`\<`T`\>

Defined in: [dimred/TSNE.js:41](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L41)

#### Parameters

##### X

`T`

The high-dimensional data.

##### parameters?

`Partial`\<[`ParametersTSNE`](../interfaces/ParametersTSNE.md)\>

Object containing parameterization of the DR method.

#### Returns

`TSNE`\<`T`\>

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

### \_gains

> **\_gains**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/TSNE.js:85](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L85)

***

### \_is\_initialized

> **\_is\_initialized**: `boolean`

Defined in: [dimred/DR.js:26](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L26)

#### Inherited from

`DR._is_initialized`

***

### \_iter

> **\_iter**: `number`

Defined in: [dimred/TSNE.js:54](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L54)

***

### \_N

> **\_N**: `number`

Defined in: [dimred/DR.js:22](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L22)

#### Inherited from

`DR._N`

***

### \_P

> **\_P**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/TSNE.js:137](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L137)

***

### \_parameters

> **\_parameters**: [`ParametersTSNE`](../interfaces/ParametersTSNE.md)

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

### \_ystep

> **\_ystep**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/TSNE.js:84](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L84)

***

### X

> **X**: [`Matrix`](Matrix.md)

Defined in: [dimred/DR.js:48](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L48)

#### Inherited from

`DR.X`

***

### Y

> **Y**: [`Matrix`](Matrix.md)

Defined in: [dimred/TSNE.js:56](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L56)

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

### check\_init()

> **check\_init**(): `DR`\<`T`, [`ParametersTSNE`](../interfaces/ParametersTSNE.md)\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersTSNE`](../interfaces/ParametersTSNE.md)\>

#### Inherited from

`DR.check_init`

***

### generator()

> **generator**(`iterations?`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/TSNE.js:157](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L157)

#### Parameters

##### iterations?

`number` = `500`

Number of iterations. Default is `500`

#### Returns

`Generator`\<`T`, `T`, `void`\>

- The projection.

#### Overrides

`DR.generator`

***

### init()

> `abstract` **init**(): `TSNE`\<`T`\>

Defined in: [dimred/TSNE.js:59](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L59)

#### Returns

`TSNE`\<`T`\>

#### Overrides

`DR.init`

***

### parameter()

#### Call Signature

> **parameter**(): [`ParametersTSNE`](../interfaces/ParametersTSNE.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersTSNE`](../interfaces/ParametersTSNE.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersTSNE`](../interfaces/ParametersTSNE.md)\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersTSNE`](../interfaces/ParametersTSNE.md)

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersTSNE`](../interfaces/ParametersTSNE.md)\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `TSNE`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersTSNE`](../interfaces/ParametersTSNE.md)

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersTSNE`](../interfaces/ParametersTSNE.md)\[`K`\]

Value of the parameter to set.

##### Returns

`TSNE`\<`T`\>

##### Inherited from

`DR.parameter`

***

### transform()

> **transform**(`iterations?`): `T`

Defined in: [dimred/TSNE.js:145](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L145)

#### Parameters

##### iterations?

`number` = `500`

Number of iterations. Default is `500`

#### Returns

`T`

The projection.

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

Defined in: [dimred/TSNE.js:270](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L270)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersTSNE`](../interfaces/ParametersTSNE.md)\>

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

***

### transform()

> `static` **transform**\<`T`\>(`X`, `parameters?`): `T`

Defined in: [dimred/TSNE.js:259](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L259)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersTSNE`](../interfaces/ParametersTSNE.md)\>

#### Returns

`T`

#### Overrides

`DR.transform`

***

### transform\_async()

> `static` **transform\_async**\<`T`\>(`X`, `parameters?`): `Promise`\<`T`\>

Defined in: [dimred/TSNE.js:282](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/TSNE.js#L282)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersTSNE`](../interfaces/ParametersTSNE.md)\>

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
