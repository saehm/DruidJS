[@saehrimnir/druidjs](../globals.md) / EigenArgs

# Interface: EigenArgs

Defined in: [linear_algebra/index.js:11](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/linear_algebra/index.js#L11)

## Properties

### max_iterations?

> `optional` **max_iterations**: `number`

Defined in: [linear_algebra/index.js:12](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/linear_algebra/index.js#L12)

The number of maxiumum iterations the algorithm should run. Default is `100`

---

### qr?

> `optional` **qr**: [`QRDecomposition`](../type-aliases/QRDecomposition.md)

Defined in: [linear_algebra/index.js:14](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/linear_algebra/index.js#L14)

The QR technique to use. Default is `qr_gramschmidt`

---

### seed?

> `optional` **seed**: `number` \| [`Randomizer`](../classes/Randomizer.md)

Defined in: [linear_algebra/index.js:13](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/linear_algebra/index.js#L13)

The seed value or a randomizer used in the algorithm. Default is `1212`

---

### tol?

> `optional` **tol**: `number`

Defined in: [linear_algebra/index.js:15](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/linear_algebra/index.js#L15)

Tolerated error for stopping criteria. Default is `1e-8`
