[@saehrimnir/druidjs](../globals.md) / SAMMON

# Class: SAMMON\<T\>

Defined in: [dimred/SAMMON.js:23](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L23)

Sammon's Mapping

A nonlinear dimensionality reduction technique that minimizes a stress
function based on the ratio of pairwise distances in high and low dimensional spaces.

## Template

## Extends

- `DR`

## Type Parameters

### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new SAMMON**\<`T`\>(`X`, `parameters?`): `SAMMON`\<`T`\>

Defined in: [dimred/SAMMON.js:35](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L35)

SAMMON's Mapping

#### Parameters

##### X

`T`

The high-dimensional data.

##### parameters?

`Partial`\<[`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>\>

Object containing parameterization of the DR
method.

#### Returns

`SAMMON`\<`T`\>

#### See

[https://arxiv.org/pdf/2009.01512.pdf](https://arxiv.org/pdf/2009.01512.pdf)

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

> **\_parameters**: [`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>

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

### distance_matrix

> **distance_matrix**: [`Matrix`](Matrix.md) \| `undefined`

Defined in: [dimred/SAMMON.js:25](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L25)

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

### \_step()

> **\_step**(): [`Matrix`](Matrix.md)

Defined in: [dimred/SAMMON.js:110](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L110)

#### Returns

[`Matrix`](Matrix.md)

---

### check_init()

> **check_init**(): `DR`\<`T`, [`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>\>

#### Inherited from

`DR.check_init`

---

### generator()

> **generator**(`max_iter?`): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/SAMMON.js:98](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L98)

Transforms the inputdata `X` to dimenionality 2.

#### Parameters

##### max_iter?

`number` = `200`

Maximum number of iteration steps. Default is `200`

#### Returns

`Generator`\<`T`, `T`, `void`\>

A generator yielding the intermediate steps of the projection of
`X`.

#### Overrides

`DR.generator`

---

### init()

> **init**(`D`): `asserts D is Matrix`

Defined in: [dimred/SAMMON.js:56](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L56)

Initializes the projection.

#### Parameters

##### D

[`Matrix`](Matrix.md) | `undefined`

#### Returns

`asserts D is Matrix`

#### Overrides

`DR.init`

---

### parameter()

#### Call Signature

> **parameter**(): [`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` _extends_ keyof [`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `SAMMON`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` _extends_ keyof [`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>\[`K`\]

Value of the parameter to set.

##### Returns

`SAMMON`\<`T`\>

##### Inherited from

`DR.parameter`

---

### transform()

> **transform**(`max_iter?`): `T`

Defined in: [dimred/SAMMON.js:82](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L82)

Transforms the inputdata `X` to dimensionality 2.

#### Parameters

##### max_iter?

`number` = `200`

Maximum number of iteration steps. Default is `200`

#### Returns

`T`

The projection of `X`.

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

Defined in: [dimred/SAMMON.js:178](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L178)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>\>

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

---

### transform()

> `static` **transform**\<`T`\>(`X`, `parameters?`): `T`

Defined in: [dimred/SAMMON.js:167](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L167)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>\>

#### Returns

`T`

#### Overrides

`DR.transform`

---

### transform_async()

> `static` **transform_async**\<`T`\>(`X`, `parameters?`): `Promise`\<`T`\>

Defined in: [dimred/SAMMON.js:190](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/SAMMON.js#L190)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters?

`Partial`\<[`ParametersSAMMON`](../interfaces/ParametersSAMMON.md)\<`AvailableInit`\>\>

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
