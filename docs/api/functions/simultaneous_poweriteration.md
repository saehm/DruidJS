[@saehrimnir/druidjs](../globals.md) / simultaneous_poweriteration

# Function: simultaneous_poweriteration()

> **simultaneous_poweriteration**(`A`, `k?`, `parameters?`): `object`

Defined in: [linear_algebra/simultaneous_poweriteration.js:19](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/linear_algebra/simultaneous_poweriteration.js#L19)

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
