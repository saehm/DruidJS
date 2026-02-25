[@saehrimnir/druidjs](../globals.md) / wasserstein

# Function: wasserstein()

> **wasserstein**(`a`, `b`): `number`

Defined in: metrics/wasserstein.js:10

Computes the 1D Wasserstein distance (Earth Mover's Distance) between two distributions.

## Parameters

### a

First distribution (histogram or probability mass)

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

### b

Second distribution (histogram or probability mass)

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

## Returns

`number`

The Wasserstein/EMD distance between `a` and `b`.

## See

[https://en.wikipedia.org/wiki/Wasserstein_metric](https://en.wikipedia.org/wiki/Wasserstein_metric)
