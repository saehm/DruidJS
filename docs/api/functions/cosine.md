[@saehrimnir/druidjs](../globals.md) / cosine

# Function: cosine()

> **cosine**(`a`, `b`): `number`

Defined in: [metrics/cosine.js:14](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/metrics/cosine.js#L14)

Computes the cosine distance (not similarity) between `a` and `b`.

## Parameters

### a

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

### b

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

## Returns

`number`

The cosine distance between `a` and `b`.

## Example

```ts
import { cosine } from "@saehrimnir/druidjs";
const a = [1, 2, 3];
const b = [4, 5, 6];
const distance = cosine(a, b); // 0.9746318461970762
```
