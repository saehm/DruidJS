<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    // IRIS embeds into a span of roughly 40, so the default threshold of 10 covers over half of
    // all point pairs and the redraw ends up repelling cluster-mates. See "Choosing
    // low_dist_thres" below -- this is a property of the algorithm, not of this port.
    const dr = new druid.LocalMAP(mistle.IRIS.values, { low_dist_thres: 25 });
    const localmap = dr.transform()
</script>

# LocalMAP

LocalMAP is a variant of [PaCMAP](/dimred/pacmap) that improves local cluster separation by dynamically resampling further pairs (FP) in the third optimization phase using nearby points in the current low-dimensional embedding, rather than keeping the random non-neighbors chosen at the start.

## How It Works

LocalMAP runs identically to PaCMAP for the first two phases, and for the first step of the third. After that it switches the nearest-neighbor gradient to a locally scaled form — attraction is multiplied by `low_dist_thres / (2 × √d_ij)`, strengthening it for pairs already close in the embedding and weakening it for far ones — and every tenth iteration it redraws the further pairs from points that lie *within* `low_dist_thres` in the current layout.

That second part is the key: a point far away in the input that has drifted close in the embedding is exactly what the repulsive term needs to push apart, and the static random set chosen at initialization rarely contains it. Rows that find no candidate inside the threshold keep the partner they had.

## Why or When to Use

Use LocalMAP when PaCMAP produces clusters that are still somewhat merged or when fine-grained local separation is important. It adds negligible overhead over PaCMAP and generally produces crisper cluster boundaries.

## Example

<RoundScatterplot :data=localmap :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const localmap = new druid.LocalMAP(data, {
    n_neighbors: 10,
    low_dist_thres: 10  // distance threshold for local FP resampling in phase 3
});

// 2. Compute the projection (450 iterations across 3 phases by default)
const projection = localmap.transform();

// Alternatively, use a generator for animation:
// for (const proj of localmap.generator()) { ... }
```

## Choosing `low_dist_thres`

`low_dist_thres` is an **absolute** distance in the embedding, and it does two jobs at once: it is
the radius within which a further pair may be redrawn, and it sets the phase 3 attraction scale to
`low_dist_thres / 2`. Both only make sense relative to how large the embedding actually is, and the
default of `10` is calibrated for the large datasets LocalMAP was designed for, whose embeddings
span far more than that.

On a small dataset the default is too small on both counts. IRIS embeds into a span of about 40, so
a threshold of 10 marks over half of all point pairs as candidates for repulsion — and since only
the `n_neighbors` nearest are excluded, most of those are cluster-mates. The clusters get pushed
apart from the inside. Measured on IRIS, mean cluster separation against the default:

| `low_dist_thres` | separation ratio |
| ---: | ---: |
| 3 | 1.5 |
| 10 (default) | 2.2 |
| 15 | 4.5 |
| 25 | 8.5 |
| PaCMAP, for comparison | 5.9 |

So LocalMAP only beats PaCMAP here once the threshold is raised past the point where it stops
targeting cluster-mates. **This is the reference algorithm's behaviour, not a quirk of this
implementation** — the numbers above are reproduced from `pacmap.LocalMAP` to within 3%.

The rule of thumb: `low_dist_thres` should be a small fraction of the embedding's span, not a
large one. Run [PaCMAP](/dimred/pacmap) first, look at the spread of the result, and scale from
there.

## Larger Datasets

LocalMAP takes the same `knn` parameter as [PaCMAP](/dimred/pacmap#larger-datasets) and defaults to
the same exact search:

```javascript
const knn = new druid.HNSW(data, { metric: druid.euclidean, ef: 100 });
const localmap = new druid.LocalMAP(data, { n_neighbors: 10, knn });
```
