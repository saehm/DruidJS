[@saehrimnir/druidjs](../globals.md) / DisjointSet

# Class: DisjointSet\<T\>

Defined in: [datastructure/DisjointSet.js:15](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/DisjointSet.js#L15)

## Template

## See

[https://en.wikipedia.org/wiki/Disjoint-set_data_structure](https://en.wikipedia.org/wiki/Disjoint-set_data_structure)

## Type Parameters

### T

`T`

## Constructors

### Constructor

> **new DisjointSet**\<`T`\>(`elements?`): `DisjointSet`\<`T`\>

Defined in: [datastructure/DisjointSet.js:19](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/DisjointSet.js#L19)

#### Parameters

##### elements?

`T`[] | `null`

#### Returns

`DisjointSet`\<`T`\>

## Methods

### find()

> **find**(`x`): `T` \| `null`

Defined in: [datastructure/DisjointSet.js:49](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/DisjointSet.js#L49)

#### Parameters

##### x

`T`

#### Returns

`T` \| `null`

---

### get_children()

> **get_children**(`x`): `Set`\<`T`\> \| `null`

Defined in: [datastructure/DisjointSet.js:98](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/DisjointSet.js#L98)

#### Parameters

##### x

`T`

#### Returns

`Set`\<`T`\> \| `null`

---

### union()

> **union**(`x`, `y`): `DisjointSet`\<`T`\>

Defined in: [datastructure/DisjointSet.js:72](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/datastructure/DisjointSet.js#L72)

#### Parameters

##### x

`T`

##### y

`T`

#### Returns

`DisjointSet`\<`T`\>
