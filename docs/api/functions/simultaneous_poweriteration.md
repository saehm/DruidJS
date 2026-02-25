[@saehrimnir/druidjs](../globals.md) / simultaneous\_poweriteration

# Function: simultaneous\_poweriteration()

> **simultaneous\_poweriteration**(`A`, `k?`, `parameters?`): `object`

Defined in: [linear\_algebra/simultaneous\_poweriteration.js:19](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/linear_algebra/simultaneous_poweriteration.js#L19)

Computes the `k` biggest Eigenvectors and Eigenvalues from Matrix `A` with the QR-Algorithm.

## Parameters

### A

[`Matrix`](../classes/Matrix.md)

The Matrix

### k?

`number` = `2`

The number of eigenvectors and eigenvalues to compute.

### parameters?

[`EigenArgs`](../interfaces/EigenArgs.md) = `{}`

Object containing parameterization of the simultanious
  poweriteration method.

## Returns

`object`

The `k` biggest eigenvectors and eigenvalues
  of Matrix `A`.

### eigenvalues

> **eigenvalues**: `Float64Array`

### eigenvectors

> **eigenvectors**: `Float64Array`\<`ArrayBufferLike`\>[]
