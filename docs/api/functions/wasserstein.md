[@saehrimnir/druidjs](../globals.md) / wasserstein

# Function: wasserstein()

> **wasserstein**(`a`, `b`): `number`

Defined in: [metrics/wasserstein.js:10](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/metrics/wasserstein.js#L10)

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

[https://en.wikipedia.org/wiki/Wasserstein\_metric](https://en.wikipedia.org/wiki/Wasserstein_metric)
