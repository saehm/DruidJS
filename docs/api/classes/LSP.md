[@saehrimnir/druidjs](../globals.md) / LSP

# Class: LSP\<T\>

Defined in: [dimred/LSP.js:23](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L23)

Least Square Projection (LSP)

A dimensionality reduction technique that uses a small set of control points
(projected with MDS) to define the projection for the rest of the data
using a Laplacian-based optimization.

## Template

## Extends

- `DR`

## Type Parameters

### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new LSP**\<`T`\>(`X`, `parameters?`): `LSP`\<`T`\>

Defined in: [dimred/LSP.js:31](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L31)

Least Squares Projection.

#### Parameters

##### X

`T`

The high-dimensional data.

##### parameters?

`Partial`\<[`ParametersLSP`](../interfaces/ParametersLSP.md)\>

Object containing parameterization of the DR method.

#### Returns

`LSP`\<`T`\>

#### See

[https://ieeexplore.ieee.org/document/4378370](https://ieeexplore.ieee.org/document/4378370)

#### Overrides

`DR.constructor`

## Properties

### \_\_input

> **\_\_input**: `T`

Defined in: [dimred/DR.js:38](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L38)

#### Inherited from

`DR.__input`

---

### \_A

> **\_A**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/LSP.js:93](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L93)

---

### \_b

> **\_b**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/LSP.js:94](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L94)

---

### \_D

> **\_D**: `number`

Defined in: [dimred/DR.js:20](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L20)

#### Inherited from

`DR._D`

---

### \_is_initialized

> **\_is_initialized**: `boolean`

Defined in: [dimred/LSP.js:49](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L49)

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

> **\_parameters**: [`ParametersLSP`](../interfaces/ParametersLSP.md)

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

> **check_init**(): `DR`\<`T`, [`ParametersLSP`](../interfaces/ParametersLSP.md)\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersLSP`](../interfaces/ParametersLSP.md)\>

#### Inherited from

`DR.check_init`

---

### generator()

> `abstract` **generator**(...`args`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/DR.js:161](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L161)

Computes the projection.

#### Parameters

##### args

...`unknown`[]

#### Returns

`Generator`\<`T`, `T`, `void`\>

The intermediate steps of the projection.

#### Inherited from

`DR.generator`

---

### init()

> **init**(): `LSP`\<`T`\>

Defined in: [dimred/LSP.js:56](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L56)

#### Returns

`LSP`\<`T`\>

#### Overrides

`DR.init`

---

### parameter()

#### Call Signature

> **parameter**(): [`ParametersLSP`](../interfaces/ParametersLSP.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersLSP`](../interfaces/ParametersLSP.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersLSP`](../interfaces/ParametersLSP.md)\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` _extends_ keyof [`ParametersLSP`](../interfaces/ParametersLSP.md)

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersLSP`](../interfaces/ParametersLSP.md)\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `LSP`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` _extends_ keyof [`ParametersLSP`](../interfaces/ParametersLSP.md)

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersLSP`](../interfaces/ParametersLSP.md)\[`K`\]

Value of the parameter to set.

##### Returns

`LSP`\<`T`\>

##### Inherited from

`DR.parameter`

---

### transform()

> **transform**(): `T`

Defined in: [dimred/LSP.js:104](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L104)

Computes the projection.

#### Returns

`T`

Returns the projection.

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

> `static` **generator**\<`T`\>(`X`, `parameters?`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/LSP.js:133](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L133)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersLSP`](../interfaces/ParametersLSP.md)\>

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

---

### transform()

> `static` **transform**\<`T`\>(`X`, `parameters?`): `T`

Defined in: [dimred/LSP.js:122](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L122)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersLSP`](../interfaces/ParametersLSP.md)\>

#### Returns

`T`

#### Overrides

`DR.transform`

---

### transform_async()

> `static` **transform_async**\<`T`\>(`X`, `parameters?`): `Promise`\<`T`\>

Defined in: [dimred/LSP.js:145](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/LSP.js#L145)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersLSP`](../interfaces/ParametersLSP.md)\>

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
