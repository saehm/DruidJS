<script setup>
import { schemeObservable10 } from "d3";
import { onMounted, ref } from "vue";
import RoundScatterplot from "./RoundScatterplot.vue";
import { runInWorker } from "./workerPool.js";

const cities = [
  // Asia
  { name: "Tokyo", lat: 35.6895, lon: 139.6917, region: "Asia" },
  { name: "Beijing", lat: 39.9042, lon: 116.4074, region: "Asia" },
  { name: "Mumbai", lat: 19.076, lon: 72.8777, region: "Asia" },
  { name: "Singapore", lat: 1.3521, lon: 103.8198, region: "Asia" },
  { name: "Seoul", lat: 37.5665, lon: 126.978, region: "Asia" },
  { name: "Bangkok", lat: 13.7563, lon: 100.5018, region: "Asia" },
  { name: "Jakarta", lat: -6.2088, lon: 106.8456, region: "Asia" },
  { name: "Manila", lat: 14.5995, lon: 120.9842, region: "Asia" },
  { name: "New Delhi", lat: 28.6139, lon: 77.209, region: "Asia" },
  { name: "Dubai", lat: 25.2048, lon: 55.2708, region: "Asia" },
  { name: "Jerusalem", lat: 31.7683, lon: 35.2137, region: "Asia" },
  { name: "Tehran", lat: 35.6892, lon: 51.389, region: "Asia" },
  { name: "Shanghai", lat: 31.2304, lon: 121.4737, region: "Asia" },
  { name: "Kolkata", lat: 22.5726, lon: 88.3639, region: "Asia" },
  { name: "Karachi", lat: 24.8607, lon: 67.0011, region: "Asia" },
  { name: "Dhaka", lat: 23.8103, lon: 90.4125, region: "Asia" },
  { name: "Guangzhou", lat: 23.1291, lon: 113.2644, region: "Asia" },
  { name: "Riyadh", lat: 24.7136, lon: 46.6753, region: "Asia" },
  { name: "Ho Chi Minh", lat: 10.8231, lon: 106.6297, region: "Asia" },
  { name: "Taipei", lat: 25.033, lon: 121.5654, region: "Asia" },
  { name: "Bengaluru", lat: 12.9716, lon: 77.5946, region: "Asia" },
  { name: "Kuala Lumpur", lat: 3.139, lon: 101.6869, region: "Asia" },
  { name: "Hanoi", lat: 21.0285, lon: 105.8542, region: "Asia" },
  { name: "Osaka", lat: 34.6937, lon: 135.5023, region: "Asia" },

  // Europe
  { name: "London", lat: 51.5074, lon: -0.1278, region: "Europe" },
  { name: "Paris", lat: 48.8566, lon: 2.3522, region: "Europe" },
  { name: "Moscow", lat: 55.7558, lon: 37.6173, region: "Europe" },
  { name: "Berlin", lat: 52.52, lon: 13.405, region: "Europe" },
  { name: "Rome", lat: 41.9028, lon: 12.4964, region: "Europe" },
  { name: "Madrid", lat: 40.4168, lon: -3.7038, region: "Europe" },
  { name: "Athens", lat: 37.9838, lon: 23.7275, region: "Europe" },
  { name: "Vienna", lat: 48.2082, lon: 16.3738, region: "Europe" },
  { name: "Amsterdam", lat: 52.3676, lon: 4.9041, region: "Europe" },
  { name: "Stockholm", lat: 59.3293, lon: 18.0686, region: "Europe" },
  { name: "Munich", lat: 48.1351, lon: 11.582, region: "Europe" },
  { name: "Milan", lat: 45.4642, lon: 9.19, region: "Europe" },
  { name: "Barcelona", lat: 41.3851, lon: 2.1734, region: "Europe" },
  { name: "Prague", lat: 50.0755, lon: 14.4378, region: "Europe" },
  { name: "Warsaw", lat: 52.2297, lon: 21.0122, region: "Europe" },
  { name: "Budapest", lat: 47.4979, lon: 19.0402, region: "Europe" },
  { name: "Lisbon", lat: 38.7223, lon: -9.1393, region: "Europe" },
  { name: "Dublin", lat: 53.3498, lon: -6.2603, region: "Europe" },
  { name: "Brussels", lat: 50.8503, lon: 4.3517, region: "Europe" },
  { name: "Copenhagen", lat: 55.6761, lon: 12.5683, region: "Europe" },
  { name: "Oslo", lat: 59.9139, lon: 10.7522, region: "Europe" },
  { name: "Helsinki", lat: 60.1695, lon: 24.9354, region: "Europe" },
  { name: "Edinburgh", lat: 55.9533, lon: -3.1883, region: "Europe" },
  { name: "Kyiv", lat: 50.4501, lon: 30.5234, region: "Europe" },
  { name: "Istanbul", lat: 41.0082, lon: 28.9784, region: "Europe" },

  // Americas
  { name: "New York", lat: 40.7128, lon: -74.006, region: "Americas" },
  { name: "Los Angeles", lat: 34.0522, lon: -118.2437, region: "Americas" },
  { name: "Chicago", lat: 41.8781, lon: -87.6298, region: "Americas" },
  { name: "Toronto", lat: 43.651, lon: -79.347, region: "Americas" },
  { name: "Vancouver", lat: 49.2827, lon: -123.1207, region: "Americas" },
  { name: "Mexico City", lat: 19.4326, lon: -99.1332, region: "Americas" },
  { name: "Rio de Janeiro", lat: -22.9068, lon: -43.1729, region: "Americas" },
  { name: "Sao Paulo", lat: -23.5505, lon: -46.6333, region: "Americas" },
  { name: "Buenos Aires", lat: -34.6037, lon: -58.3816, region: "Americas" },
  { name: "Santiago", lat: -33.4489, lon: -70.6693, region: "Americas" },
  { name: "Bogota", lat: 4.711, lon: -74.0721, region: "Americas" },
  { name: "Lima", lat: -12.0464, lon: -77.0428, region: "Americas" },
  { name: "Montreal", lat: 45.5017, lon: -73.5673, region: "Americas" },
  { name: "Houston", lat: 29.7604, lon: -95.3698, region: "Americas" },
  { name: "Miami", lat: 25.7617, lon: -80.1918, region: "Americas" },
  { name: "Atlanta", lat: 33.749, lon: -84.388, region: "Americas" },
  { name: "Boston", lat: 42.3601, lon: -71.0589, region: "Americas" },
  { name: "San Francisco", lat: 37.7749, lon: -122.4194, region: "Americas" },
  { name: "Seattle", lat: 47.6062, lon: -122.3321, region: "Americas" },
  { name: "Washington", lat: 38.9072, lon: -77.0369, region: "Americas" },
  { name: "Dallas", lat: 32.7767, lon: -96.797, region: "Americas" },
  { name: "Caracas", lat: 10.4806, lon: -66.9036, region: "Americas" },
  { name: "Havana", lat: 23.1136, lon: -82.3666, region: "Americas" },
  { name: "Belo Horizonte", lat: -19.9167, lon: -43.9345, region: "Americas" },

  // Africa
  { name: "Cape Town", lat: -33.9249, lon: 18.4241, region: "Africa" },
  { name: "Cairo", lat: 30.0444, lon: 31.2357, region: "Africa" },
  { name: "Lagos", lat: 6.5244, lon: 3.3792, region: "Africa" },
  { name: "Nairobi", lat: -1.2921, lon: 36.8219, region: "Africa" },
  { name: "Johannesburg", lat: -26.2041, lon: 28.0473, region: "Africa" },
  { name: "Casablanca", lat: 33.5731, lon: -7.5898, region: "Africa" },
  { name: "Accra", lat: 5.6037, lon: -0.187, region: "Africa" },
  { name: "Algiers", lat: 36.7538, lon: 3.0588, region: "Africa" },
  { name: "Kinshasa", lat: -4.4419, lon: 15.2663, region: "Africa" },
  { name: "Luanda", lat: -8.8147, lon: 13.2302, region: "Africa" },
  { name: "Dar es Salaam", lat: -6.7924, lon: 39.2083, region: "Africa" },
  { name: "Khartoum", lat: 15.5007, lon: 32.5599, region: "Africa" },
  { name: "Abidjan", lat: 5.3599, lon: -4.0083, region: "Africa" },
  { name: "Alexandria", lat: 31.2001, lon: 29.9187, region: "Africa" },
  { name: "Addis Ababa", lat: 9.0227, lon: 38.7468, region: "Africa" },
  { name: "Dakar", lat: 14.7167, lon: -17.4677, region: "Africa" },
  { name: "Durban", lat: -29.8587, lon: 31.0218, region: "Africa" },
  { name: "Tunis", lat: 36.8065, lon: 10.1815, region: "Africa" },

  // Oceania
  { name: "Sydney", lat: -33.8688, lon: 151.2093, region: "Oceania" },
  { name: "Melbourne", lat: -37.8136, lon: 144.9631, region: "Oceania" },
  { name: "Brisbane", lat: -27.4698, lon: 153.0251, region: "Oceania" },
  { name: "Perth", lat: -31.9505, lon: 115.8605, region: "Oceania" },
  { name: "Auckland", lat: -36.8485, lon: 174.7633, region: "Oceania" },
  { name: "Wellington", lat: -41.2865, lon: 174.7762, region: "Oceania" },
  { name: "Adelaide", lat: -34.9285, lon: 138.6007, region: "Oceania" },
  { name: "Gold Coast", lat: -28.0167, lon: 153.4, region: "Oceania" },
  { name: "Hobart", lat: -42.8821, lon: 147.3272, region: "Oceania" },
  { name: "Darwin", lat: -12.4634, lon: 130.8456, region: "Oceania" },
  { name: "Christchurch", lat: -43.532, lon: 172.6362, region: "Oceania" },
  { name: "Suva", lat: -18.1404, lon: 178.439, region: "Oceania" },
  { name: "Port Moresby", lat: -9.4431, lon: 147.1803, region: "Oceania" },
  { name: "Noumea", lat: -22.2758, lon: 166.458, region: "Oceania" },
  { name: "Papeete", lat: -17.535, lon: -149.5696, region: "Oceania" },
];

const results = ref(null);
const loading = ref(true);

const regions = Array.from(new Set(cities.map((c) => c.region))).sort();

function getRegionColor(region) {
  const index = regions.indexOf(region);
  return schemeObservable10[index % schemeObservable10.length];
}

async function run() {
  loading.value = true;
  const toRad = (d) => (d * Math.PI) / 180;
  const data = cities.map((c) => [toRad(c.lat), toRad(c.lon)]);

  try {
    const projection = await runInWorker(
      "DR",
      "SMACOF",
      data,
      {
        metric: "haversine",
      },
      500,
    );
    results.value = projection;
  } catch (e) {
    console.error("Error in Geospatial UMAP:", e);
  } finally {
    loading.value = false;
  }
}

onMounted(() => {
  run();
});
</script>

<template>
  <div class="geospatial-projection">
    <div class="header">
      <p class="description">
        Projecting world cities using the spherical Haversine metric and SMACOF.
      </p>
      <button @click="run" :disabled="loading" class="btn">Regenerate</button>
    </div>

    <div class="content">
      <div class="plot-box">
        <div v-if="loading" class="loading-overlay">
          <div class="spinner"></div>
          Projecting...
        </div>
        <RoundScatterplot
          v-if="results"
          :data="results"
          :labels="cities.map((c) => c.name)"
          :colors="cities.map((c) => getRegionColor(c.region))"
          :size="400"
          :radius="8"
          :style="{ opacity: loading ? 0.3 : 1, transition: 'opacity 0.3s' }"
        />
      </div>

      <div class="sidebar">
        <div class="legend">
          <h3>Regions</h3>
          <div class="legend-items">
            <div v-for="region in regions" :key="region" class="legend-item">
              <div class="color-swatch" :style="{ backgroundColor: getRegionColor(region) }"></div>
              <span>{{ region }}</span>
            </div>
          </div>
        </div>

        <div class="info">
          <p>
            This projection uses the <b>Haversine Distance</b> to correctly calculate distances
            along the Earth's curved surface.
          </p>
          <p><em>Hover over points in the scatterplot to see the city names!</em></p>
        </div>
      </div>
    </div>
  </div>
</template>

<style scoped>
.geospatial-projection {
  border: 1px solid var(--vp-c-divider);
  border-radius: 12px;
  padding: 24px;
  background: var(--vp-c-bg-soft);
  margin: 20px 0;
}

.header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 20px;
  flex-wrap: wrap;
  gap: 10px;
}

.description {
  margin: 0;
  color: var(--vp-c-text-2);
  font-size: 0.95em;
}

.btn {
  padding: 6px 16px;
  border-radius: 6px;
  border: 1px solid var(--vp-c-divider);
  background: var(--vp-c-bg);
  cursor: pointer;
  font-weight: 500;
  transition: border-color 0.2s;
}

.btn:hover:not(:disabled) {
  border-color: var(--vp-c-brand);
}

.btn:disabled {
  opacity: 0.6;
  cursor: not-allowed;
}

.content {
  display: flex;
  flex-wrap: wrap;
  gap: 1em;
  align-items: stretch;
}

.plot-box {
  position: relative;
  min-height: 400px;
  min-width: 400px;
  display: flex;
  justify-content: center;
  align-items: center;
  border-radius: 8px;
  background: var(--vp-c-bg);
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
}

.loading-overlay {
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 12px;
  color: var(--vp-c-text-2);
  z-index: 10;
}

.spinner {
  width: 30px;
  height: 30px;
  border: 3px solid var(--vp-c-divider);
  border-top-color: var(--vp-c-brand);
  border-radius: 50%;
  animation: spin 1s linear infinite;
}

.sidebar {
  flex: 1;
  min-width: 150px;
  display: flex;
  flex-direction: column;
  gap: 20px;
}

.legend {
  background: var(--vp-c-bg);
  padding: 16px;
  border-radius: 8px;
  border: 1px solid var(--vp-c-divider);
}

.legend h3 {
  margin: 0 0 16px 0;
  font-size: 1.1em;
  font-weight: 600;
}

.legend-items {
  display: flex;
  flex-direction: column;
  gap: 12px;
}

.legend-item {
  display: flex;
  align-items: center;
  gap: 12px;
  font-size: 0.9em;
  color: var(--vp-c-text-1);
}

.color-swatch {
  width: 14px;
  height: 14px;
  border-radius: 50%;
  flex-shrink: 0;
}

.info {
  background: var(--vp-c-bg);
  padding: 16px;
  border-radius: 8px;
  border: 1px solid var(--vp-c-divider);
  font-size: 0.9em;
  color: var(--vp-c-text-2);
}

.info p {
  margin: 0 0 10px 0;
}

.info p:last-child {
  margin-bottom: 0;
}

@keyframes spin {
  to {
    transform: rotate(360deg);
  }
}
</style>
