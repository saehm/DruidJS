[@saehrimnir/druidjs](../globals.md) / TopoMap

# Class: TopoMap\<T\>

Defined in: [dimred/TopoMap.js:21](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L21)

TopoMap

A 0-dimensional Homology Preserving Projection of High-Dimensional Data.
It aims to preserve the topological structure of the data by maintaining
the connectivity of a minimum spanning tree.

## Template

## Extends

- `DR`

## Type Parameters

### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

## Constructors

### Constructor

> **new TopoMap**\<`T`\>(`X`, `parameters`): `TopoMap`\<`T`\>

Defined in: [dimred/TopoMap.js:29](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L29)

TopoMap: A 0-dimensional Homology Preserving Projection of High-Dimensional Data.

#### Parameters

##### X

`T`

The high-dimensional data.

##### parameters

`Partial`\<[`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)\>

Object containing parameterization of the DR method.

#### Returns

`TopoMap`\<`T`\>

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

### \_disjoint_set

> **\_disjoint_set**: [`DisjointSet`](DisjointSet.md)\<`Float64Array`\<`ArrayBufferLike`\>\> \| `undefined`

Defined in: [dimred/TopoMap.js:66](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L66)

---

### \_distance_matrix

> **\_distance_matrix**: [`Matrix`](Matrix.md)

Defined in: [dimred/TopoMap.js:32](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L32)

---

### \_Emst

> **\_Emst**: `number`[][] \| `undefined`

Defined in: [dimred/TopoMap.js:94](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L94)

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

> **\_parameters**: [`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)

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

> **check_init**(): `DR`\<`T`, [`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)\>

Defined in: [dimred/DR.js:202](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L202)

If the respective DR method has an `init` function, call it before `transform`.

#### Returns

`DR`\<`T`, [`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)\>

#### Inherited from

`DR.check_init`

---

### generator()

> **generator**(): `Generator`\<`T`, `T`, `void`\>

Defined in: [dimred/TopoMap.js:322](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L322)

Transforms the inputdata `X` to dimensionality 2.

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

---

### init()

> **init**(): `TopoMap`\<`T`\>

Defined in: [dimred/TopoMap.js:91](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L91)

Initializes TopoMap. Sets all projcted points to zero, and computes a minimum spanning tree.

#### Returns

`TopoMap`\<`T`\>

#### Overrides

`DR.init`

---

### parameter()

#### Call Signature

> **parameter**(): [`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)

Defined in: [dimred/DR.js:74](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L74)

Get all Parameters.

##### Returns

[`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`): [`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)\[`K`\]

Defined in: [dimred/DR.js:80](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L80)

Get value of given parameter.

##### Type Parameters

###### K

`K` _extends_ keyof [`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)

##### Parameters

###### name

`K`

Name of the parameter.

##### Returns

[`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)\[`K`\]

##### Inherited from

`DR.parameter`

#### Call Signature

> **parameter**\<`K`\>(`name`, `value`): `TopoMap`\<`T`\>

Defined in: [dimred/DR.js:87](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/DR.js#L87)

Set value of given parameter.

##### Type Parameters

###### K

`K` _extends_ keyof [`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)

##### Parameters

###### name

`K`

Name of the parameter.

###### value

[`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)\[`K`\]

Value of the parameter to set.

##### Returns

`TopoMap`\<`T`\>

##### Inherited from

`DR.parameter`

---

### transform()

> **transform**(): `T`

Defined in: [dimred/TopoMap.js:290](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L290)

Transforms the inputdata `X` to dimensionality 2.

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

Defined in: [dimred/TopoMap.js:366](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L366)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters

`Partial`\<[`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)\>

#### Returns

`Generator`\<`T`, `T`, `void`\>

#### Overrides

`DR.generator`

---

### transform()

> `static` **transform**\<`T`\>(`X`, `parameters`): `T`

Defined in: [dimred/TopoMap.js:355](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L355)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters

`Partial`\<[`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)\>

#### Returns

`T`

#### Overrides

`DR.transform`

---

### transform_async()

> `static` **transform_async**\<`T`\>(`X`, `parameters`): `Promise`\<`T`\>

Defined in: [dimred/TopoMap.js:378](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/TopoMap.js#L378)

#### Type Parameters

##### T

`T` _extends_ [`InputType`](../type-aliases/InputType.md)

#### Parameters

##### X

`T`

##### parameters

`Partial`\<[`ParametersTopoMap`](../interfaces/ParametersTopoMap.md)\>

#### Returns

`Promise`\<`T`\>

#### Overrides

`DR.transform_async`
