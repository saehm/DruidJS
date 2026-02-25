import * as druid from "../../dist/druid.js";

// Store persistent instances to avoid rebuilding large indices (e.g. KNN)
const instances = new Map();

self.onmessage = async (e) => {
  const { id, task, method, data, params = {}, iterations, instanceId } = e.data;
  let result;

  try {
    const Algorithm = druid[method];
    if (task !== "Search" && !Algorithm) throw new Error(`Algorithm ${method} not found`);

    if (task !== "Search") {
      if (!params.metric) {
        params.metric = druid.euclidean;
      } else if (typeof params.metric === "string") {
        params.metric = druid[params.metric];
      }
    }

    if (task === "DR" || task === "AnimatedDR") {
      const dr = new Algorithm(data, params);
      if (dr.init) dr.init();

      if (task === "DR") {
        result = dr.transform(iterations);
        self.postMessage({ id, result, success: true });
      } else {
        const generator = dr.generator(iterations || 200);
        let step = 0;
        for (const projection of generator) {
          self.postMessage({ id, result: projection, step, success: true, done: false });
          step++;
        }
        self.postMessage({ id, success: true, done: true });
      }
    } else if (task === "Clustering") {
      // KMeans specifically takes { K: k }
      const clusteringParams =
        method === "KMeans" || method === "KMedoids" ? { K: params.k || 3 } : params;
      const clustering = new Algorithm(data, clusteringParams);
      if (clustering.init) clustering.init();
      result = clustering.get_cluster_list(params.cut_value);
      self.postMessage({ id, result, success: true });
    } else if (task === "BuildKNN") {
      const knn = new Algorithm(data, params);
      if (instanceId) {
        instances.set(instanceId, knn);
        if (instances.size > 5) {
          const firstKey = instances.keys().next().value;
          instances.delete(firstKey);
        }
      }
      self.postMessage({ id, success: true });
    } else if (task === "Search") {
      const knn = instanceId ? instances.get(instanceId) : null;
      if (!knn) throw new Error(`Instance ${instanceId} not found`);
      const { query, k = 10 } = params;
      result = knn.search(query, k);
      self.postMessage({ id, result, success: true });
    }
  } catch (error) {
    self.postMessage({
      id,
      error: `Worker Error [${task}:${method}]: ${error.message}`,
      success: false,
    });
  }
};
