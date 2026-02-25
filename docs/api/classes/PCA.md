[@saehrimnir/druidjs](../globals.md) / PCA

# Class: PCA\<T\>

Defined in: [dimred/PCA.js:29](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L29)

Principal Component Analysis (PCA)

A linear dimensionality reduction technique that identifies the axes (principal components)
along which the variance of the data is maximized.

## Template

## See

[MDS](MDS.md) for another linear alternative

## Example

```ts
import * as druid from "@saehrimnir/druidjs";

const X = [[1, 2], [3, 4], [5, 6]];
const pca = new druid.PCA(X, { d: 2 });
const Y = pca.transform();
// [[x1, y1], [x2, y2], [x3, y3]]
```

## Extends

- `DR`

## Type Parameters

### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new PCA**\<`T`\>(`X`, `parameters?`): `PCA`\<`T`\>

Defined in: [dimred/PCA.js:34](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L34)

#### Parameters

##### X

`T`

The high-dimensional data.

##### parameters?

`Partial`\<[`ParametersPCA`](../interfaces/ParametersPCA.md)\> = `{}`

Object containing parameterization of the DR method.

#### Returns

`PCA`\<`T`\>

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

> **\_parameters**: [`ParametersPCA`](../interfaces/ParametersPCA.md)

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

### V

> **V**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/PCA.js:90](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L90)

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

### check\_init()

> **check\_init**(): `DR`\<`T`, [`ParametersPCA`](../interfaces/ParametersPCA.md)\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersPCA`](../interfaces/ParametersPCA.md)\>

#### Inherited from

`DR.check_init`

***

### generator()

> **generator**(): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/PCA.js:47](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L47)

Transforms the inputdata `X` to dimensionality `d`.

#### Returns

`Generator`\<`T`, `T`, `void`\>

A generator yielding the intermediate steps of the projection.

#### Overrides

`DR.generator`

***

### init()

> `abstract` **init**(...`args`): `void`

Defined in: [dimred/DR.js:193](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L193)

#### Parameters

##### args

...`unknown`[]

#### Returns

`void`

#### Inherited from

`DR.init`

***

### parameter()

#### Call Signature

> **parameter**(): [`ParametersPCA`](../interfaces/ParametersPCA.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersPCA`](../interfaces/ParametersPCA.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersPCA`](../interfaces/ParametersPCA.md)\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersPCA`](../interfaces/ParametersPCA.md)

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersPCA`](../interfaces/ParametersPCA.md)\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `PCA`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersPCA`](../interfaces/ParametersPCA.md)

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersPCA`](../interfaces/ParametersPCA.md)\[`K`\]

Value of the parameter to set.

##### Returns

`PCA`\<`T`\>

##### Inherited from

`DR.parameter`

***

### principal\_components()

> **principal\_components**(): [`Matrix`](Matrix.md)

Defined in: [dimred/PCA.js:80](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L80)

Computes the `d` principal components of Matrix `X`.

#### Returns

[`Matrix`](Matrix.md)

***

### transform()

> **transform**(): `T`

Defined in: [dimred/PCA.js:57](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L57)

Transforms the inputdata `X` to dimensionality `d`.

#### Returns

`T`

- The projected data.

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

Defined in: [dimred/PCA.js:111](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L111)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersPCA`](../interfaces/ParametersPCA.md)\>

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

***

### principal\_components()

> `static` **principal\_components**\<`T`\>(`X`, `parameters`): [`Matrix`](Matrix.md)

Defined in: [dimred/PCA.js:100](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L100)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters

`Partial`\<[`ParametersPCA`](../interfaces/ParametersPCA.md)\>

#### Returns

[`Matrix`](Matrix.md)

***

### transform()

> `static` **transform**\<`T`\>(`X`, `parameters`): `T`

Defined in: [dimred/PCA.js:70](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L70)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters

`Partial`\<[`ParametersPCA`](../interfaces/ParametersPCA.md)\>

#### Returns

`T`

#### Overrides

`DR.transform`

***

### transform\_async()

> `static` **transform\_async**\<`T`\>(`X`, `parameters?`): `Promise`\<`T`\>

Defined in: [dimred/PCA.js:123](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/PCA.js#L123)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersPCA`](../interfaces/ParametersPCA.md)\>

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
