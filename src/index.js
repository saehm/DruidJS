import pkg from "../package.json" with { type: "json" };

const version = pkg.version;
export { version };

/** @import {Matrix} from "./matrix/index.js" */
/** @typedef {Matrix | Float64Array[] | number[][]} InputType*/

//export { version } from "../package.json" with { type: "json" };
export * from "./clustering/index.js";
export * from "./datastructure/index.js";
export * from "./dimred/index.js";
export * from "./knn/index.js";
export * from "./linear_algebra/index.js";
export * from "./matrix/index.js";
export * from "./metrics/index.js";
export * from "./numerical/index.js";
export * from "./optimization/index.js";
export * from "./util/index.js";
