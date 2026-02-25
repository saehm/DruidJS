[@saehrimnir/druidjs](../globals.md) / ParametersOptics

# Interface: ParametersOptics

Defined in: [clustering/index.js:35](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L35)

## Properties

### epsilon

> **epsilon**: `number`

Defined in: [clustering/index.js:36](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L36)

The minimum distance which defines whether a point is a neighbor or not.

***

### metric

> **metric**: [`Metric`](../type-aliases/Metric.md)

Defined in: [clustering/index.js:38](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L38)

The distance metric which defines the distance between two points of the points. Default is `euclidean`

***

### min\_points

> **min\_points**: `number`

Defined in: [clustering/index.js:37](https://github.com/saehm/DruidJS/blob/dcefb767556f3c4cf43bdeb36d2750ae8d5ef93c/src/clustering/index.js#L37)

The minimum number of points which a point needs to create a cluster. (Should be higher than 1, else each point creates a cluster.)
