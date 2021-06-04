import resolve from "@rollup/plugin-node-resolve";
import { terser } from "rollup-plugin-terser";
import jsdoc from 'rollup-plugin-jsdoc';
import json from '@rollup/plugin-json';
import { babel } from '@rollup/plugin-babel';
import meta from "./package.json";


const copyright = `// ${meta.homepage} v${meta.version} Copyright ${(new Date).getFullYear()} ${meta.author.name}`;

export default [
  {
    input: "index.js",
    output: {
      extend: true,
      banner: copyright,
      file: "dist/druid.js",
      format: "umd",
      indent: false,
      name: "druid"
    },
    plugins: [
        json({
            compact: true,
            exclude: 'node_modules/**',
        }),
        resolve(),
        jsdoc({
            args: ["-r", "-d", "docs"],
            config: "jsdoc.config.json",
        }),
    ]
  },
  {
    input: "index.js",
    plugins: [
        json({
            compact: true
        }),
        resolve(),
        terser({
            format: {
                preamble: copyright
            }
        })
    ],
    output: {
      extend: true,
      file: "dist/druid.min.js",
      format: "umd",
      indent: false,
      name: "druid"
    }
  },
  {
    input: "index.js",
    plugins: [
        json({
            compact: true
        }),
        resolve(),
        terser({
            format: {
                comments: "all",
                preamble: copyright
            },
            keep_classnames: true,
            keep_fnames: true,
            compress: true,
        })
    ],
    output: {
      extend: true,
      file: "dist/druid.esm.js",
      format: "es",
      indent: false,
      name: "druid"
    }
  },
  {
    input: "index.js",
    plugins: [
        json({
            compact: true
        }),
        babel({ 
            babelHelpers: 'runtime',
            plugins: [
                ["@babel/plugin-proposal-nullish-coalescing-operator"],
                ["@babel/plugin-transform-runtime", {loose: true}]
            ]
        }),
        resolve(),
        terser({
            format: {
                preamble: copyright
            },
            keep_classnames: true,
            keep_fnames: true,
            compress: true,
        })
    ],
    external: [/@babel\/runtime/],
    output: {
      extend: true,
      file: "dist/druid.es6.js",
      format: "umd",
      indent: false,
      name: "druid"
    }
  }
];