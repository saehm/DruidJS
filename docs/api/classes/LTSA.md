[@saehrimnir/druidjs](../globals.md) / LTSA

# Class: LTSA\<T\>

Defined in: [dimred/LTSA.js:22](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LTSA.js#L22)

Local Tangent Space Alignment (LTSA)

A nonlinear dimensionality reduction algorithm that represents the local
geometry of the manifold by tangent spaces and then aligns them to reveal
the global structure.

## Template

## Extends

- `DR`

## Type Parameters

### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new LTSA**\<`T`\>(`X`, `parameters`): `LTSA`\<`T`\>

Defined in: [dimred/LTSA.js:30](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LTSA.js#L30)

Local Tangent Space Alignment

#### Parameters

##### X

`T`

The high-dimensional data.

##### parameters

`Partial`\<[`ParametersLTSA`](../interfaces/ParametersLTSA.md)\>

Object containing parameterization of the DR method.

#### Returns

`LTSA`\<`T`\>

#### See

[https://epubs.siam.org/doi/abs/10.1137/S1064827502419154](https://epubs.siam.org/doi/abs/10.1137/S1064827502419154)

#### Overrides

`DR.constructor`

## Properties

### \_\_input

> **\_\_input**: `T`

Defined in: [dimred/DR.js:38](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L38)

#### Inherited from

`DR.__input`

---

### \_D

> **\_D**: `number`

Defined in: [dimred/DR.js:20](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L20)

#### Inherited from

`DR._D`

---

### \_is_initialized

> **\_is_initialized**: `boolean`

Defined in: [dimred/DR.js:26](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L26)

#### Inherited from

`DR._is_initialized`

---

### \_N

> **\_N**: `number`

Defined in: [dimred/DR.js:22](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L22)

#### Inherited from

`DR._N`

---

### \_parameters

> **\_parameters**: [`ParametersLTSA`](../interfaces/ParametersLTSA.md)

Defined in: [dimred/DR.js:41](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L41)

#### Inherited from

`DR._parameters`

---

### \_randomizer

> **\_randomizer**: [`Randomizer`](Randomizer.md)

Defined in: [dimred/DR.js:24](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L24)

#### Inherited from

`DR._randomizer`

---

### \_type

> **\_type**: `"array"` \| `"matrix"` \| `"typed"`

Defined in: [dimred/DR.js:46](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L46)

#### Inherited from

`DR._type`

---

### X

> **X**: [`Matrix`](Matrix.md)

Defined in: [dimred/DR.js:48](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L48)

#### Inherited from

`DR.X`

---

### Y

> **Y**: [`Matrix`](Matrix.md)

Defined in: [dimred/DR.js:50](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L50)

#### Inherited from

`DR.Y`

## Accessors

### projection

#### Get Signature

> **get** **projection**(): `T`

Defined in: [dimred/DR.js:211](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L211)

##### Returns

`T`

The projection in the type of input `X`.

#### Inherited from

`DR.projection`

## Methods

### check_init()

> **check_init**(): `DR`\<`T`, [`ParametersLTSA`](../interfaces/ParametersLTSA.md)\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersLTSA`](../interfaces/ParametersLTSA.md)\>

#### Inherited from

`DR.check_init`

---

### generator()

> **generator**(): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/LTSA.js:63](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LTSA.js#L63)

Transforms the inputdata `X` to dimensionality `d`.

#### Returns

`Generator`\<`T`, `T`, `void`\>

A generator yielding the intermediate steps of the projection.

#### Overrides

`DR.generator`

---

### init()

> `abstract` **init**(...`args`): `void`

Defined in: [dimred/DR.js:193](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L193)

#### Parameters

##### args

...`unknown`[]

#### Returns

`void`

#### Inherited from

`DR.init`

---

### parameter()

#### Call Signature

> **parameter**(): [`ParametersLTSA`](../interfaces/ParametersLTSA.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersLTSA`](../interfaces/ParametersLTSA.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersLTSA`](../interfaces/ParametersLTSA.md)\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` _extends_ keyof [`ParametersLTSA`](../interfaces/ParametersLTSA.md)

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersLTSA`](../interfaces/ParametersLTSA.md)\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `LTSA`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` _extends_ keyof [`ParametersLTSA`](../interfaces/ParametersLTSA.md)

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersLTSA`](../interfaces/ParametersLTSA.md)\[`K`\]

Value of the parameter to set.

##### Returns

`LTSA`\<`T`\>

##### Inherited from

`DR.parameter`

---

### transform()

> **transform**(): `T`

Defined in: [dimred/LTSA.js:73](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LTSA.js#L73)

Transforms the inputdata `X` to dimenionality `d`.

#### Returns

`T`

#### Overrides

`DR.transform`

---

### transform_async()

> **transform_async**(...`args`): `Promise`\<`T`\>

Defined in: [dimred/DR.js:233](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L233)

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

---

### generator()

> `static` **generator**\<`T`\>(`X`, `parameters`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/LTSA.js:131](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LTSA.js#L131)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters

`Partial`\<[`ParametersLTSA`](../interfaces/ParametersLTSA.md)\>

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

---

### transform()

> `static` **transform**\<`T`\>(`X`, `parameters`): `T`

Defined in: [dimred/LTSA.js:120](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LTSA.js#L120)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters

`Partial`\<[`ParametersLTSA`](../interfaces/ParametersLTSA.md)\>

#### Returns

`T`

#### Overrides

`DR.transform`

---

### transform_async()

> `static` **transform_async**\<`T`\>(`X`, `parameters`): `Promise`\<`T`\>

Defined in: [dimred/LTSA.js:143](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LTSA.js#L143)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters

`Partial`\<[`ParametersLTSA`](../interfaces/ParametersLTSA.md)\>

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
