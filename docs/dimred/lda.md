<script setup>
    import RoundScatterplot from "../components/RoundScatterplot.vue";

    import * as druid from "../../dist/druid.js";
    import * as mistle from "@saehrimnir/mistle";
    
    const labels = mistle.IRIS.labels;
    const dr = new druid.LDA(mistle.IRIS.values, {labels: mistle.IRIS.labels});
    const lda = dr.transform()
</script>

# LDA

Linear Discriminant Analysis (LDA) is a technique used to find a linear combination of features that characterizes or separates two or more classes of objects or events.

## How It Works

Linear Discriminant Analysis (LDA) finds a linear combination of features that maximizes the ratio of between-class variance to within-class variance in the dataset.

## Why or When to Use

Use LDA for supervised dimensionality reduction where class labels are known, aiming to separate distinct classes as much as possible, often as a preprocessing step for classification.

## Example

<RoundScatterplot :data=lda :labels="labels" />

## How-to (Code)

```javascript
import * as druid from "@saehrimnir/druidjs";

const data = [
  /* ... multi-dimensional data ... */
];
const classLabels = [
  /* ... labels ... */
];

// 1. Initialize the algorithm
const lda = new druid.LDA(data, { labels: classLabels });

// 2. Compute the projection
const projection = lda.transform();
```
