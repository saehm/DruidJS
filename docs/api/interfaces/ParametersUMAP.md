[@saehrimnir/druidjs](../globals.md) / ParametersUMAP

# Interface: ParametersUMAP

Defined in: [dimred/index.js:149](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L149)

## Properties

### \_initial_alpha?

> `optional` **\_initial_alpha**: `number`

Defined in: [dimred/index.js:160](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L160)

The initial learning rate for the optimization.

---

### \_n_epochs?

> `optional` **\_n_epochs**: `number`

Defined in: [dimred/index.js:159](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L159)

The number of training epochs.

---

### \_negative_sample_rate?

> `optional` **\_negative_sample_rate**: `number`

Defined in: [dimred/index.js:158](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L158)

The number of negative samples per positive sample.

---

### \_repulsion_strength?

> `optional` **\_repulsion_strength**: `number`

Defined in: [dimred/index.js:157](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L157)

Weighting applied to negative samples.

---

### \_set_op_mix_ratio?

> `optional` **\_set_op_mix_ratio**: `number`

Defined in: [dimred/index.js:156](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L156)

Interpolate between union and intersection.

---

### \_spread?

> `optional` **\_spread**: `number`

Defined in: [dimred/index.js:155](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L155)

The effective scale of embedded points.

---

### d?

> `optional` **d**: `number`

Defined in: [dimred/index.js:153](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L153)

the dimensionality of the projection.

---

### local_connectivity?

> `optional` **local_connectivity**: `number`

Defined in: [dimred/index.js:151](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L151)

number of nearest neighbors connected in the local neighborhood.

---

### metric?

> `optional` **metric**: [`Metric`](../type-aliases/Metric.md) \| `"precomputed"`

Defined in: [dimred/index.js:154](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L154)

the metric which defines the distance between two points in the high-dimensional space.

---

### min_dist?

> `optional` **min_dist**: `number`

Defined in: [dimred/index.js:152](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L152)

controls how tightly points get packed together.

---

### n_neighbors?

> `optional` **n_neighbors**: `number`

Defined in: [dimred/index.js:150](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L150)

size of the local neighborhood.

---

### seed?

> `optional` **seed**: `number`

Defined in: [dimred/index.js:161](https://github.com/saehm/DruidJS/blob/a8c3d973d427068e4ee01a4685c87a4d0f8b20c6/src/dimred/index.js#L161)

the seed for the random number generator.
