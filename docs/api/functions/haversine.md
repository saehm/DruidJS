[@saehrimnir/druidjs](../globals.md) / haversine

# Function: haversine()

> **haversine**(`a`, `b`): `number`

Defined in: [metrics/haversine.js:10](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/metrics/haversine.js#L10)

Computes the Haversine distance between two points on a sphere of unit length 1. Multiply the result with the radius of the sphere. (For instance Earth's radius is 6371km)

## Parameters

### a

Point [lat1, lon1] in radians

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

### b

Point [lat2, lon2] in radians

`number`[] | `Float64Array`\<`ArrayBufferLike`\>

## Returns

`number`

The Haversine distance between `a` and `b`.

## See

[https://en.wikipedia.org/wiki/Haversine\_formula](https://en.wikipedia.org/wiki/Haversine_formula)
