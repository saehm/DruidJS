[@saehrimnir/druidjs](../globals.md) / SQDMDS

# Class: SQDMDS\<T\>

Defined in: [dimred/SQDMDS.js:21](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L21)

SQuadMDS (Stochastic Quartet MDS)

A lean Stochastic Quartet MDS improving global structure preservation in
neighbor embedding like t-SNE and UMAP.

## Template

## Extends

- `DR`

## Type Parameters

### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new SQDMDS**\<`T`\>(`X`, `parameters?`): `SQDMDS`\<`T`\>

Defined in: [dimred/SQDMDS.js:30](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L30)

SQuadMDS: a lean Stochastic Quartet MDS improving global structure preservation in neighbor embedding like t-SNE
and UMAP.

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)\>

#### Returns

`SQDMDS`\<`T`\>

#### See

[https://arxiv.org/pdf/2202.12087.pdf](https://arxiv.org/pdf/2202.12087.pdf)

#### Overrides

`DR.constructor`

## Properties

### \_\_input

> **\_\_input**: `T`

Defined in: [dimred/DR.js:38](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L38)

#### Inherited from

`DR.__input`

***

### \_add

> **\_add**: (...`summands`) => `Float64Array`\<`ArrayBufferLike`\> \| `undefined`

Defined in: [dimred/SQDMDS.js:54](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L54)

***

### \_D

> **\_D**: `number`

Defined in: [dimred/DR.js:20](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L20)

#### Inherited from

`DR._D`

***

### \_decay\_start

> **\_decay\_start**: `number` \| `undefined`

Defined in: [dimred/SQDMDS.js:96](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L96)

***

### \_distance\_exaggeration

> **\_distance\_exaggeration**: `boolean` \| `undefined`

Defined in: [dimred/SQDMDS.js:137](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L137)

***

### \_grads

> **\_grads**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/SQDMDS.js:63](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L63)

***

### \_HD\_metric

> **\_HD\_metric**: (`i`, `j`, `X`) => `number` \| `undefined`

Defined in: [dimred/SQDMDS.js:73](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L73)

***

### \_HD\_metric\_exaggeration

> **\_HD\_metric\_exaggeration**: (`i`, `j`, `X`) => `number` \| `undefined`

Defined in: [dimred/SQDMDS.js:75](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L75)

***

### \_indices

> **\_indices**: `number`[] \| `undefined`

Defined in: [dimred/SQDMDS.js:64](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L64)

***

### \_is\_initialized

> **\_is\_initialized**: `boolean`

Defined in: [dimred/DR.js:26](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L26)

#### Inherited from

`DR._is_initialized`

***

### \_LR

> **\_LR**: `number` \| `undefined`

Defined in: [dimred/SQDMDS.js:59](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L59)

***

### \_LR\_init

> **\_LR\_init**: `number` \| `undefined`

Defined in: [dimred/SQDMDS.js:58](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L58)

***

### \_minus

> **\_minus**: (`a`, `b`) => `Float64Array`\<`ArrayBufferLike`\> \| `undefined`

Defined in: [dimred/SQDMDS.js:56](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L56)

***

### \_momentums

> **\_momentums**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/SQDMDS.js:62](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L62)

***

### \_mult

> **\_mult**: (`a`, `v`) => `Float64Array`\<`ArrayBufferLike`\> \| `undefined`

Defined in: [dimred/SQDMDS.js:57](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L57)

***

### \_N

> **\_N**: `number`

Defined in: [dimred/DR.js:22](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L22)

#### Inherited from

`DR._N`

***

### \_offset

> **\_offset**: `number` \| `undefined`

Defined in: [dimred/SQDMDS.js:61](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L61)

***

### \_parameters

> **\_parameters**: [`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)

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

### \_sub\_div

> **\_sub\_div**: (`x`, `y`, `div`) => `Float64Array`\<`ArrayBufferLike`\> \| `undefined`

Defined in: [dimred/SQDMDS.js:55](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L55)

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

### \_\_add()

> **\_\_add**(`d`): (...`summands`) => `Float64Array`\<`ArrayBufferLike`\>

Defined in: [dimred/SQDMDS.js:425](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L425)

Inline!

#### Parameters

##### d

`number`

#### Returns

> (...`summands`): `Float64Array`\<`ArrayBufferLike`\>

##### Parameters

###### summands

...`Float64Array`\<`ArrayBufferLike`\>[]

##### Returns

`Float64Array`\<`ArrayBufferLike`\>

***

### \_\_minus()

> **\_\_minus**(`d`): (`a`, `b`) => `Float64Array`\<`ArrayBufferLike`\>

Defined in: [dimred/SQDMDS.js:411](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L411)

Inline!

#### Parameters

##### d

`number`

#### Returns

> (`a`, `b`): `Float64Array`\<`ArrayBufferLike`\>

##### Parameters

###### a

`Float64Array`\<`ArrayBufferLike`\>

###### b

`Float64Array`\<`ArrayBufferLike`\>

##### Returns

`Float64Array`\<`ArrayBufferLike`\>

***

### \_\_mult()

> **\_\_mult**(`d`): (`a`, `v`) => `Float64Array`\<`ArrayBufferLike`\>

Defined in: [dimred/SQDMDS.js:444](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L444)

Inline!

#### Parameters

##### d

`number`

#### Returns

> (`a`, `v`): `Float64Array`\<`ArrayBufferLike`\>

##### Parameters

###### a

`Float64Array`\<`ArrayBufferLike`\>

###### v

`number`

##### Returns

`Float64Array`\<`ArrayBufferLike`\>

***

### \_\_sub\_div()

> **\_\_sub\_div**(`d`): (`x`, `y`, `div`) => `Float64Array`\<`ArrayBufferLike`\>

Defined in: [dimred/SQDMDS.js:458](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L458)

Creates a new array `(x - y) / div`.

#### Parameters

##### d

`number`

#### Returns

> (`x`, `y`, `div`): `Float64Array`\<`ArrayBufferLike`\>

##### Parameters

###### x

`Float64Array`\<`ArrayBufferLike`\>

###### y

`Float64Array`\<`ArrayBufferLike`\>

###### div

`number`

##### Returns

`Float64Array`\<`ArrayBufferLike`\>

***

### \_fill\_MDS\_grads()

> **\_fill\_MDS\_grads**(`Y`, `grads`, `exaggeration?`, `zero_grad?`): [`Matrix`](Matrix.md)

Defined in: [dimred/SQDMDS.js:205](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L205)

Computes the gradients.

#### Parameters

##### Y

[`Matrix`](Matrix.md)

The Projection.

##### grads

[`Matrix`](Matrix.md)

The gradients.

##### exaggeration?

`boolean` = `false`

Whether or not to use early exaggeration. Default is `false`

##### zero\_grad?

`boolean` = `true`

Whether or not to reset the gradient in the beginning. Default is `true`

#### Returns

[`Matrix`](Matrix.md)

The gradients.

***

### check\_init()

> **check\_init**(): `DR`\<`T`, [`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)\>

#### Inherited from

`DR.check_init`

***

### generator()

> **generator**(`iterations?`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/SQDMDS.js:109](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L109)

Computes the projection.

#### Parameters

##### iterations?

`number` = `500`

Number of iterations. Default is `500`

#### Returns

`Generator`\<`T`, `T`, `void`\>

The intermediate steps of the projection.

#### Overrides

`DR.generator`

***

### init()

> `abstract` **init**(): `void`

Defined in: [dimred/SQDMDS.js:49](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L49)

#### Returns

`void`

#### Overrides

`DR.init`

***

### parameter()

#### Call Signature

> **parameter**(): [`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `SQDMDS`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` *extends* keyof [`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)\[`K`\]

Value of the parameter to set.

##### Returns

`SQDMDS`\<`T`\>

##### Inherited from

`DR.parameter`

***

### transform()

> **transform**(`iterations?`): `T`

Defined in: [dimred/SQDMDS.js:93](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L93)

Computes the projection.

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

Defined in: [dimred/SQDMDS.js:481](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L481)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)\>

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

***

### transform()

> `static` **transform**\<`T`\>(`X`, `parameters?`): `T`

Defined in: [dimred/SQDMDS.js:470](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L470)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)\>

#### Returns

`T`

#### Overrides

`DR.transform`

***

### transform\_async()

> `static` **transform\_async**\<`T`\>(`X`, `parameters?`): `Promise`\<`T`\>

Defined in: [dimred/SQDMDS.js:493](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/dimred/SQDMDS.js#L493)

#### Type Parameters

##### T

`T` *extends* [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersSQDMDS`](../interfaces/ParametersSQDMDS.md)\>

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
