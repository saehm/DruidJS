[@saehrimnir/druidjs](../globals.md) / LDA

# Class: LDA\<T\>

Defined in: [dimred/LDA.js:20](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/LDA.js#L20)

Linear Discriminant Analysis (LDA)

A supervised dimensionality reduction technique that finds the axes that
maximize the separation between multiple classes.

## Template

## Extends

- `DR`

## Type Parameters

### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new LDA**\<`T`\>(`X`, `parameters`): `LDA`\<`T`\>

Defined in: [dimred/LDA.js:28](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/LDA.js#L28)

Linear Discriminant Analysis.

#### Parameters

##### X

`T`

The high-dimensional data.

##### parameters

`Partial`\<[`ParametersLDA`](../interfaces/ParametersLDA.md)\> & `object`

Object containing parameterization of the DR method.

#### Returns

`LDA`\<`T`\>

#### See

[https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x](https://onlinelibrary.wiley.com/doi/10.1111/j.1469-1809.1936.tb02137.x)

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

> **\_parameters**: [`ParametersLDA`](../interfaces/ParametersLDA.md)

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

> **check\_init**(): `DR`\<`T`, [`ParametersLDA`](../interfaces/ParametersLDA.md)\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersLDA`](../interfaces/ParametersLDA.md)\>

#### Inherited from

`DR.check_init`

***

### generator()

> **generator**(): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/LDA.js:41](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/LDA.js#L41)

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

> **parameter**(): [`ParametersLDA`](../interfaces/ParametersLDA.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersLDA`](../interfaces/ParametersLDA.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersLDA`](../interfaces/ParametersLDA.md)\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersLDA`](../interfaces/ParametersLDA.md)

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersLDA`](../interfaces/ParametersLDA.md)\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `LDA`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersLDA`](../interfaces/ParametersLDA.md)

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersLDA`](../interfaces/ParametersLDA.md)\[`K`\]

Value of the parameter to set.

##### Returns

`LDA`\<`T`\>

##### Inherited from

`DR.parameter`

***

### transform()

> **transform**(): `T`

Defined in: [dimred/LDA.js:51](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/LDA.js#L51)

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

> `static` **generator**\<`T`, `Para`\>(`X`, `parameters`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/LDA.js:137](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/LDA.js#L137)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

##### Para

`Para` *extends* `object`

#### Parameters

##### X

`T`

##### parameters

`Para`

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

***

### transform()

> `static` **transform**\<`T`, `Para`\>(`X`, `parameters`): `T`

Defined in: [dimred/LDA.js:124](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/LDA.js#L124)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

##### Para

`Para` *extends* `object`

#### Parameters

##### X

`T`

##### parameters

`Para`

#### Returns

`T`

#### Overrides

`DR.transform`

***

### transform\_async()

> `static` **transform\_async**\<`T`, `Para`\>(`X`, `parameters`): `Promise`\<`T`\>

Defined in: [dimred/LDA.js:151](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/LDA.js#L151)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

##### Para

`Para` *extends* `object`

#### Parameters

##### X

`T`

##### parameters

`Para`

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
