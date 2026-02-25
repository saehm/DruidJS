export class DruidWorker {
  constructor() {
    this.worker = new Worker(new URL("./druid.worker.js", import.meta.url), { type: "module" });
    this.callbacks = new Map();
    this.streams = new Map();
    this.nextId = 0;

    this.worker.onmessage = (e) => {
      const { id, result, error, success, done, step } = e.data;

      const streamCallback = this.streams.get(id);
      if (streamCallback) {
        if (!success) {
          streamCallback.reject(new Error(error));
          this.streams.delete(id);
        } else if (done) {
          streamCallback.resolve();
          this.streams.delete(id);
        } else {
          streamCallback.onProgress(result, step);
        }
        return;
      }

      const callback = this.callbacks.get(id);
      if (callback) {
        if (success) {
          callback.resolve(result);
        } else {
          callback.reject(new Error(error));
        }
        this.callbacks.delete(id);
      }
    };
  }

  run(task, method, data, params = {}, iterations = undefined, instanceId = undefined) {
    const id = this.nextId++;
    return new Promise((resolve, reject) => {
      this.callbacks.set(id, { resolve, reject });
      this.worker.postMessage({ id, task, method, data, params, iterations, instanceId });
    });
  }

  runStream(
    task,
    method,
    data,
    params = {},
    iterations = undefined,
    onProgress,
    instanceId = undefined,
  ) {
    const id = this.nextId++;
    return new Promise((resolve, reject) => {
      this.streams.set(id, { resolve, reject, onProgress });
      this.worker.postMessage({ id, task, method, data, params, iterations, instanceId });
    });
  }

  terminate() {
    this.worker.terminate();
  }
}

const pool = [];
const MAX_WORKERS = navigator.hardwareConcurrency || 4;

function getWorker() {
  let druidWorker = pool.find((w) => w.callbacks.size === 0 && w.streams.size === 0);
  if (!druidWorker && pool.length < MAX_WORKERS) {
    druidWorker = new DruidWorker();
    pool.push(druidWorker);
  }
  return druidWorker || pool[Math.floor(Math.random() * pool.length)];
}

export async function runInWorker(task, method, data, params = {}, iterations = undefined) {
  return getWorker().run(task, method, data, params, iterations);
}

export async function runStreamInWorker(
  task,
  method,
  data,
  params = {},
  iterations = undefined,
  onProgress,
) {
  return getWorker().runStream(task, method, data, params, iterations, onProgress);
}

export async function createPersistentInstance(method, data, params = {}) {
  const instanceId = Math.random().toString(36).substring(7);
  const worker = getWorker();
  await worker.run("BuildKNN", method, data, params, undefined, instanceId);
  return { instanceId, worker };
}

export async function usePersistentInstance(instance, task, params = {}) {
  return instance.worker.run(task, undefined, undefined, params, undefined, instance.instanceId);
}
