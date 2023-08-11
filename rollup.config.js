import resolve from "@rollup/plugin-node-resolve";
import terser from "@rollup/plugin-terser";
import json from "@rollup/plugin-json";

const onwarn = (message, warn) => {
    if (message.code === "CIRCULAR_DEPENDENCY") return;
    warn(message);
};

export default [
    {
        input: "src/index.js",
        output: {
            sourcemap: true,
            extend: true,
            file: "dist/druid.js",
            format: "umd",
            indent: false,
            name: "druid",
        },
        plugins: [
            json({
                compact: true,
                exclude: "node_modules/**",
            }),
            resolve(),
        ],
        onwarn,
    },
    {
        input: "src/index.js",
        plugins: [
            json({
                compact: true,
            }),
            resolve(),
            terser(),
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
    } /* 
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
    }, */,
];
