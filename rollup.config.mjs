import json from "@rollup/plugin-json";
import resolve from "@rollup/plugin-node-resolve";
import { dts } from "rollup-plugin-dts";

export default [
    {
        input: "src/index.js",
        output: [
            {
                file: "dist/druid.js",
                format: "esm",
                sourcemap: true,
            },
            {
                file: "dist/druid.cjs",
                format: "cjs",
                sourcemap: true,
            },
        ],
        plugins: [resolve(), json()],
    },
    {
        input: "dist/types/index.d.ts",
        output: [
            {
                file: "dist/druid.d.ts",
                format: "es",
                sourcemap: true,
            },
            {
                file: "dist/druid.d.cts",
                format: "es",
                sourcemap: true,
            },
        ],
        plugins: [dts()],
    },
];
