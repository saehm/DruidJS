[@saehrimnir/druidjs](../globals.md) / Randomizer

# Class: Randomizer

Defined in: [util/randomizer.js:7](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L7)

## Constructors

### Constructor

> **new Randomizer**(`_seed?`): `Randomizer`

Defined in: [util/randomizer.js:28](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L28)

Mersenne Twister random number generator.

#### Parameters

##### \_seed?

`number`

The seed for the random number generator. If `_seed == null` then
the actual time gets used as seed. Default is `new Date().getTime()`

#### Returns

`Randomizer`

#### See

https://github.com/bmurray7/mersenne-twister-examples/blob/master/javascript-mersenne-twister.js

## Properties

### \_LOWER_MASK

> **\_LOWER_MASK**: `number` = `0x7fffffff`

Defined in: [util/randomizer.js:12](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L12)

---

### \_M

> **\_M**: `number` = `397`

Defined in: [util/randomizer.js:9](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L9)

---

### \_MATRIX_A

> **\_MATRIX_A**: `number` = `0x9908b0df`

Defined in: [util/randomizer.js:10](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L10)

---

### \_mt

> **\_mt**: `number`[]

Defined in: [util/randomizer.js:15](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L15)

---

### \_mti

> **\_mti**: `number`

Defined in: [util/randomizer.js:17](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L17)

---

### \_N

> **\_N**: `number` = `624`

Defined in: [util/randomizer.js:8](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L8)

---

### \_seed

> **\_seed**: `number`

Defined in: [util/randomizer.js:19](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L19)

---

### \_UPPER_MASK

> **\_UPPER_MASK**: `number` = `0x80000000`

Defined in: [util/randomizer.js:11](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L11)

---

### \_val

> **\_val**: `number` \| `null` \| `undefined`

Defined in: [util/randomizer.js:113](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L113)

## Accessors

### random

#### Get Signature

> **get** **random**(): `number`

Defined in: [util/randomizer.js:63](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L63)

Returns a float between 0 and 1.

##### Returns

`number`

- A random number between [0, 1]

---

### random_int

#### Get Signature

> **get** **random_int**(): `number`

Defined in: [util/randomizer.js:72](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L72)

Returns an integer between 0 and MAX_INTEGER.

##### Returns

`number`

- A random integer.

---

### seed

#### Get Signature

> **get** **seed**(): `number`

Defined in: [util/randomizer.js:54](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L54)

Returns the seed of the random number generator.

##### Returns

`number`

- The seed.

#### Set Signature

> **set** **seed**(`_seed`): `void`

Defined in: [util/randomizer.js:36](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L36)

##### Parameters

###### \_seed

`number`

##### Returns

`void`

## Methods

### choice()

> **choice**\<`T`\>(`A`, `n`): `T`[]

Defined in: [util/randomizer.js:132](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L132)

#### Type Parameters

##### T

`T`

Returns samples from an input Matrix or Array.

#### Parameters

##### A

`T`[]

The input Matrix or Array.

##### n

`number`

The number of samples.

#### Returns

`T`[]

A random selection form `A` of `n` samples.

---

### gauss_random()

> **gauss_random**(): `number`

Defined in: [util/randomizer.js:109](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L109)

#### Returns

`number`

---

### choice()

> `static` **choice**\<`T`\>(`A`, `n`, `seed?`): `T`[]

Defined in: [util/randomizer.js:171](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/util/randomizer.js#L171)

#### Type Parameters

##### T

`T`

Returns samples from an input Matrix or Array.

#### Parameters

##### A

`T`[]

The input Matrix or Array.

##### n

`number`

The number of samples.

##### seed?

`number` = `1212`

The seed for the random number generator.

#### Returns

`T`[]

- A random selection form `A` of `n` samples.
