# Global Spatial Projections

This example demonstrates how dimensionality reduction can be used for **geospatial analysis** by using a specialized metric.

Most dimensionality reduction techniques use `euclidean` distance by default, which works fine for flat planes. However, when working with coordinates on the Earth's surface, the Curvature of the Earth becomes significant.

Using the **Haversine Distance** (Great-Circle Distance) ensures that the distance between two points is calculated correctly along the sphere.

<script setup>
import GeospatialProjection from '../components/GeospatialProjection.vue'
</script>

<GeospatialProjection />

## Why use Haversine?

- **Spherical Accuracy**: Traditional Euclidean distance "tunnels" through the Earth, while Haversine follows the surface.
- **Metric-Aware DR**: DruidJS allows you to swap metrics easily. Here, we use [`SMACOF`](/api/classes/SMACOF) (Metric MDS) with `metric: 'haversine'`.
- **Topological Clusters**: Notice how cities cluster by continent (Oceania, Europe, Americas, etc.) naturally based on their real-world proximity.

## City Dataset

We are projecting 100 major world cities based solely on their latitude and longitude coordinates.

## How-to (Code)

To use the `haversine` metric, ensure your coordinates are in **radians**. You can then pass the `haversine` metric function to algorithms like [`UMAP`](/api/classes/UMAP) or [`MDS`](/api/classes/MDS).

```javascript
import * as druid from "@saehrimnir/druidjs";

// Helper to convert degrees to radians
const toRad = (d) => (d * Math.PI) / 180;

// Cities [lat, lon] in radians
const data = [
  [toRad(35.68), toRad(139.69)], // Tokyo
  [toRad(40.71), toRad(-74.01)], // New York
  // ... more cities
];

// Initialize SMACOF with the haversine metric
const smacof = new druid.SMACOF(data, {
  metric: druid.haversine,
  d: 2,
});

// Run the iterative optimization
const projection = smacof.transform(500);
```
