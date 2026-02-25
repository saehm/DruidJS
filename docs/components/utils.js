export function applyZScore(data) {
  if (!data || data.length === 0) return data;
  const n = data.length;
  const m = data[0].length || 0;
  if (m === 0) return data;

  // Deep copy
  const result = new Array(n);
  for (let i = 0; i < n; i++) {
    result[i] = new Array(m);
  }

  for (let j = 0; j < m; j++) {
    let mean = 0;
    for (let i = 0; i < n; i++) {
      mean += data[i][j];
    }
    mean /= n;

    let std = 0;
    for (let i = 0; i < n; i++) {
      std += (data[i][j] - mean) ** 2;
    }
    std = Math.sqrt(std / n);

    if (std === 0) std = 1; // Prevent division by zero

    for (let i = 0; i < n; i++) {
      result[i][j] = (data[i][j] - mean) / std;
    }
  }
  return result;
}
