export function max(values) {
  let max;
  for (const value of values) {
    if (value != null
        && (max < value || (max === undefined && value >= value))) {
      max = value;
    }
  }
  return max;
}
