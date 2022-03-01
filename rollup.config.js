import resolve from "@rollup/plugin-node-resolve";
import { terser } from "rollup-plugin-terser";
import json from "@rollup/plugin-json";
import meta from "./package.json";
import ts from "rollup-plugin-ts";

const copyright = `// ${meta.homepage} v${meta.version} Copyright ${new Date().getFullYear()} ${meta.author.name}`;

const onwarn = (message, warn) => {
    if (message.code === "CIRCULAR_DEPENDENCY") return;
    warn(message);
};

export default [
    {
        input: "index.js",
        output: {
            sourcemap: true,
            extend: true,
            banner: copyright,
            file: "dist/druid.js",
            format: "umd",
            indent: false,
            name: "druid",
        },
        plugins: [
            ts(),
            json({
                compact: true,
                exclude: "node_modules/**",
            }),
            resolve()
        ],
        onwarn,
    },
    {
        input: "index.js",
        plugins: [
            json({
                compact: true,
            }),
            resolve(),
            terser({
                format: {
                    preamble: copyright,
                },
            }),
        ],
        output: {
            sourcemap: true,
            extend: true,
            file: "dist/druid.min.js",
            format: "umd",
            indent: false,
            name: "druid",
        },
        onwarn,
    },
    {
        input: "index.js",
        plugins: [
            json({
                compact: true,
            }),
            resolve(),
            terser({
                format: {
                    comments: "all",
                    preamble: copyright,
                },
                keep_classnames: true,
                keep_fnames: true,
                compress: true,
            }),
        ],
        output: {
            sourcemap: true,
            extend: true,
            file: "dist/druid.esm.js",
            format: "es",
            indent: false,
            name: "druid",
        },
        onwarn,
    },
];
