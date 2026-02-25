[@saehrimnir/druidjs](../globals.md) / UMAP

# Class: UMAP\<T\>

Defined in: [dimred/UMAP.js:41](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L41)

Uniform Manifold Approximation and Projection (UMAP)

A novel manifold learning technique for dimensionality reduction. UMAP is constructed
from a theoretical framework based on Riemannian geometry and algebraic topology.
It is often faster than t-SNE while preserving more of the global structure.

## Template

## See

 - [Paper](https://arxiv.org/abs/1802.03426|UMAP)
 - [TSNE](TSNE.md) for a similar visualization technique

## Example

```ts
import * as druid from "@saehrimnir/druidjs";

const X = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]];
const umap = new druid.UMAP(X, {
    n_neighbors: 15,
    min_dist: 0.1,
    d: 2,
    seed: 42
});

const Y = umap.transform(500); // 500 iterations
// [[x1, y1], [x2, y2], [x3, y3]]
```

## Extends

- `DR`

## Type Parameters

### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new UMAP**\<`T`\>(`X`, `parameters?`): `UMAP`\<`T`\>

Defined in: [dimred/UMAP.js:46](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L46)

#### Parameters

##### X

`T`

The high-dimensional data.

##### parameters?

`Partial`\<[`ParametersUMAP`](../interfaces/ParametersUMAP.md)\>

Object containing parameterization of the DR method.

#### Returns

`UMAP`\<`T`\>

#### Overrides

`DR.constructor`

## Properties

### \_\_input

> **\_\_input**: `T`

Defined in: [dimred/DR.js:38](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L38)

#### Inherited from

`DR.__input`

***

### \_a

> **\_a**: `number` \| `undefined`

Defined in: [dimred/UMAP.js:331](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L331)

***

### \_alpha

> **\_alpha**: `number` \| `undefined`

Defined in: [dimred/UMAP.js:479](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L479)

***

### \_b

> **\_b**: `number` \| `undefined`

Defined in: [dimred/UMAP.js:332](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L332)

***

### \_D

> **\_D**: `number`

Defined in: [dimred/DR.js:20](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L20)

#### Inherited from

`DR._D`

***

### \_epoch\_of\_next\_negative\_sample

> **\_epoch\_of\_next\_negative\_sample**: `Float32Array`\<`ArrayBuffer`\> \| `undefined`

Defined in: [dimred/UMAP.js:341](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L341)

***

### \_epoch\_of\_next\_sample

> **\_epoch\_of\_next\_sample**: `Float32Array`\<`ArrayBuffer`\> \| `undefined`

Defined in: [dimred/UMAP.js:340](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L340)

***

### \_epochs\_per\_negative\_sample

> **\_epochs\_per\_negative\_sample**: `Float32Array`\<`ArrayBuffer`\> \| `undefined`

Defined in: [dimred/UMAP.js:339](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L339)

***

### \_epochs\_per\_sample

> **\_epochs\_per\_sample**: `Float32Array`\<`ArrayBufferLike`\> \| `undefined`

Defined in: [dimred/UMAP.js:338](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L338)

***

### \_graph

> **\_graph**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/UMAP.js:333](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L333)

***

### \_head

> **\_head**: `number`[] \| `undefined`

Defined in: [dimred/UMAP.js:335](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L335)

***

### \_is\_initialized

> **\_is\_initialized**: `boolean`

Defined in: [dimred/DR.js:26](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L26)

#### Inherited from

`DR._is_initialized`

***

### \_iter

> **\_iter**: `number`

Defined in: [dimred/UMAP.js:82](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L82)

***

### \_N

> **\_N**: `number`

Defined in: [dimred/DR.js:22](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L22)

#### Inherited from

`DR._N`

***

### \_parameters

> **\_parameters**: [`ParametersUMAP`](../interfaces/ParametersUMAP.md)

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

### \_tail

> **\_tail**: `number`[] \| `undefined`

Defined in: [dimred/UMAP.js:336](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L336)

***

### \_type

> **\_type**: `"array"` \| `"matrix"` \| `"typed"`

Defined in: [dimred/DR.js:46](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L46)

#### Inherited from

`DR._type`

***

### \_weights

> **\_weights**: `number`[] \| `undefined`

Defined in: [dimred/UMAP.js:337](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L337)

***

### X

> **X**: [`Matrix`](Matrix.md)

Defined in: [dimred/DR.js:48](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L48)

#### Inherited from

`DR.X`

***

### Y

> **Y**: [`Matrix`](Matrix.md)

Defined in: [dimred/UMAP.js:84](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L84)

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

> **check\_init**(): `DR`\<`T`, [`ParametersUMAP`](../interfaces/ParametersUMAP.md)\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersUMAP`](../interfaces/ParametersUMAP.md)\>

#### Inherited from

`DR.check_init`

***

### generator()

> **generator**(`iterations?`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/UMAP.js:370](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L370)

#### Parameters

##### iterations?

`number` = `350`

Number of iterations. Default is `350`

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

***

### graph()

> **graph**(): `object`

Defined in: [dimred/UMAP.js:345](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L345)

#### Returns

`object`

##### cols

> **cols**: `number`[] \| `undefined`

##### rows

> **rows**: `number`[] \| `undefined`

##### weights

> **weights**: `number`[] \| `undefined`

***

### init()

> **init**(): `UMAP`\<`T`\>

Defined in: [dimred/UMAP.js:324](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L324)

Computes all necessary

#### Returns

`UMAP`\<`T`\>

#### Overrides

`DR.init`

***

### parameter()

#### Call Signature

> **parameter**(): [`ParametersUMAP`](../interfaces/ParametersUMAP.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersUMAP`](../interfaces/ParametersUMAP.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersUMAP`](../interfaces/ParametersUMAP.md)\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersUMAP`](../interfaces/ParametersUMAP.md)

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersUMAP`](../interfaces/ParametersUMAP.md)\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `UMAP`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersUMAP`](../interfaces/ParametersUMAP.md)

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersUMAP`](../interfaces/ParametersUMAP.md)\[`K`\]

Value of the parameter to set.

##### Returns

`UMAP`\<`T`\>

##### Inherited from

`DR.parameter`

***

### transform()

> **transform**(`iterations?`): `T`

Defined in: [dimred/UMAP.js:354](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L354)

#### Parameters

##### iterations?

`number` = `350`

Number of iterations. Default is `350`

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

Defined in: [dimred/UMAP.js:502](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L502)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersUMAP`](../interfaces/ParametersUMAP.md)\>

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

***

### transform()

> `static` **transform**\<`T`\>(`X`, `parameters?`): `T`

Defined in: [dimred/UMAP.js:491](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L491)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersUMAP`](../interfaces/ParametersUMAP.md)\>

#### Returns

`T`

#### Overrides

`DR.transform`

***

### transform\_async()

> `static` **transform\_async**\<`T`\>(`X`, `parameters?`): `Promise`\<`T`\>

Defined in: [dimred/UMAP.js:514](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/UMAP.js#L514)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersUMAP`](../interfaces/ParametersUMAP.md)\>

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
