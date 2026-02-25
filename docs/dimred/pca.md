<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.PCA(mistle.IRIS.values);
    const pca = dr.transform()
</script>

# PCA

Principal Component Analysis, or **PCA**, is a powerful dimensionality reduction technique used to simplify complex datasets while losing as little information as possible. Think of it as finding the "best angle" to take a photo of a 3D object so that you can still tell exactly what it is in a 2D picture.

## How It Works

At its core, PCA transforms a large set of variables into a smaller one that still contains most of the original information. It does this by identifying **Principal Components**:

1. **Standardization:** The data is scaled so that each variable contributes equally to the analysis.
2. **Covariance Matrix Computation:** It looks for correlations between variables to see how they vary from the mean in relation to each other.
3. **Eigenvectors and Eigenvalues:** PCA calculates new axes (Principal Components). The first principal component is the direction that captures the maximum variance (spread) in the data.
4. **Feature Vector:** It discards the components that account for very little variance, effectively "dropping" the noise and keeping the signal.

## Why Use It?

- **Data Compression:** It reduces the number of variables, making models run faster and requiring less memory.
- **Visualization:** It can crush high-dimensional data (like 100 different metrics) down to 2D or 3D so humans can actually plot and understand it.
- **Noise Reduction:** By focusing only on the directions with the most variance, it often filters out random fluctuations that don't matter.
- **Interpretability:** It provides a clear, interpretable summary of the data's structure, making it easier to understand and communicate insights.

## Example

<RoundScatterplot :data=pca :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];

// 1. Initialize the algorithm
const pca = new druid.PCA(data);

// 2. Compute the projection
const projection = pca.transform();
```
