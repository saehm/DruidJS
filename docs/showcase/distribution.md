# Earth Mover Analysis

This example demonstrates using the **Wasserstein Distance** (also known as Earth Mover's Distance) to compare probability distributions.

Unlike Euclidean distance, which compares values bin-by-bin, Wasserstein distance measures the minimum "work" required to transform one distribution into another. This makes it ideal for comparing histograms where peaks might be shifted.

In the visualization below, we generate several synthetic 1D distributions (histograms) and use **Multi-Dimensional Scaling ([MDS](/api/classes/MDS))** with the `wasserstein` metric to project them into 2D space.

<script setup>
import DistributionMDS from '../components/DistributionMDS.vue'
</script>

<DistributionMDS />

## How it works

1. **Distribution Generation**: We create diverse 1D distributions (Gaussian, Uniform, Bimodal, etc.).
2. **Wasserstein Metric**: We calculate the pairwise 1D Wasserstein distance between all histograms.
3. **MDS Projection**: MDS takes the distance matrix and finds a 2D configuration that preserves these "earth mover" distances.

::: warning
Notice how the `Uniform`, `Bimodal`, and `Wide Gaussian` distributions cluster closely together—changing one into another requires moving dirt only short, local distances. In contrast, moving the peak of `Gaussian 1` (left) all the way into `Gaussian 3` (right) requires carrying 100% of the distribution across the entire axis, which is reflected via a massive projected distance, pushing them physically far apart!
:::

## How-to (Code)

To perform Earth Mover Analysis in DruidJS, you can use the `MDS` algorithm and specify `wasserstein` as the metric. This is particularly useful for comparing histograms or probability density functions.

```javascript
import * as druid from "@saehrimnir/druidjs";

// Your data: An array of histograms (normalized distributions)
const data = [
  [0.1, 0.2, 0.4, 0.2, 0.1], // Distribution A
  [0.0, 0.1, 0.3, 0.4, 0.2], // Distribution B
  // ... more distributions
];

// Initialize MDS with the wasserstein metric
const mds = new druid.MDS(data, {
  metric: druid.wasserstein,
  d: 2,
});

// Project the distributions into 2D space
const projection = mds.transform();
```
