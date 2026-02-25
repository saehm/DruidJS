[@saehrimnir/druidjs](../globals.md) / Heap

# Class: Heap\<T\>

Defined in: [datastructure/Heap.js:8](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L8)

## Template

## Type Parameters

### T

`T`

## Constructors

### Constructor

> **new Heap**\<`T`\>(`elements?`, `accessor`, `comparator?`): `Heap`\<`T`\>

Defined in: [datastructure/Heap.js:26](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L26)

A heap is a datastructure holding its elements in a specific way, so that the top element would be the first
entry of an ordered list.

#### Parameters

##### elements?

Contains the elements for the Heap. `elements` can be null.

`T`[] | `null`

##### accessor

(`d`) => `number`

Function returns the value of the element.

##### comparator?

Function returning true or false
defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a
Max_heap). Default is `"min"`

`"max"` | `"min"` | [`Comparator`](../type-aliases/Comparator.md)

#### Returns

`Heap`\<`T`\>

#### See

[https://en.wikipedia.org/wiki/Binary_heap](https://en.wikipedia.org/wiki/Binary_heap)

## Properties

### \_accessor()

> **\_accessor**: (`d`) => `number`

Defined in: [datastructure/Heap.js:28](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L28)

#### Parameters

##### d

`T`

#### Returns

`number`

---

### \_comparator

> **\_comparator**: [`Comparator`](../type-aliases/Comparator.md)

Defined in: [datastructure/Heap.js:13](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L13)

---

### \_container

> **\_container**: `object`[]

Defined in: [datastructure/Heap.js:10](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L10)

#### element

> **element**: `T`

#### value

> **value**: `number`

## Accessors

### empty

#### Get Signature

> **get** **empty**(): `boolean`

Defined in: [datastructure/Heap.js:227](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L227)

Returns false if the the heap has entries, true if the heap has no entries.

##### Returns

`boolean`

---

### first

#### Get Signature

> **get** **first**(): \{ `element`: `T`; `value`: `number`; \} \| `null`

Defined in: [datastructure/Heap.js:171](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L171)

Returns the top entry of the heap without removing it.

##### Returns

\{ `element`: `T`; `value`: `number`; \} \| `null`

Object consists of the element and its value (computed by
`accessor`).

---

### length

#### Get Signature

> **get** **length**(): `number`

Defined in: [datastructure/Heap.js:218](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L218)

The size of the heap.

##### Returns

`number`

## Methods

### data()

> **data**(): `T`[]

Defined in: [datastructure/Heap.js:200](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L200)

Returns elements of container array.

#### Returns

`T`[]

Array consisting the elements.

---

### iterate()

> **iterate**(): `Generator`\<`T`, `void`, `unknown`\>

Defined in: [datastructure/Heap.js:180](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L180)

Yields the raw data

#### Returns

`Generator`\<`T`, `void`, `unknown`\>

#### Yields

Object consists of the element and its value (computed by `accessor`}).

---

### pop()

> **pop**(): \{ `element`: `T`; `value`: `number`; \} \| `null`

Defined in: [datastructure/Heap.js:150](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L150)

Removes and returns the top entry of the heap.

#### Returns

\{ `element`: `T`; `value`: `number`; \} \| `null`

Object consists of the element and its value (computed by
`accessor`}).

---

### push()

> **push**(`element`): `Heap`\<`T`\>

Defined in: [datastructure/Heap.js:111](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L111)

Pushes the element to the heap.

#### Parameters

##### element

`T`

#### Returns

`Heap`\<`T`\>

---

### raw_data()

> **raw_data**(): `object`[]

Defined in: [datastructure/Heap.js:209](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L209)

Returns the container array.

#### Returns

`object`[]

The container array.

---

### toArray()

> **toArray**(): `T`[]

Defined in: [datastructure/Heap.js:191](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L191)

Returns the heap as ordered array.

#### Returns

`T`[]

Array consisting the elements ordered by `comparator`.

---

### heapify()

> `static` **heapify**\<`T`\>(`elements`, `accessor`, `comparator?`): `Heap`\<`T`\>

Defined in: [datastructure/Heap.js:62](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/Heap.js#L62)

Creates a Heap from an Array

#### Type Parameters

##### T

`T`

#### Parameters

##### elements

`T`[]

Contains the elements for the Heap.

##### accessor

(`d`) => `number`

Function returns the value of the element.

##### comparator?

Function returning true or false
defining the wished order of the Heap, or String for predefined function. ("min" for a Min-Heap, "max" for a
Max_heap). Default is `"min"`

`"max"` | `"min"` | [`Comparator`](../type-aliases/Comparator.md)

#### Returns

`Heap`\<`T`\>
