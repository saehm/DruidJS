import resolve from "@rollup/plugin-node-resolve";
import { terser } from "rollup-plugin-terser";
import jsdoc from 'rollup-plugin-jsdoc';
import json from '@rollup/plugin-json';
import meta from "./package.json";

//const copyright = `// ${meta.homepage} v${meta.version} Copyright ${(new Date).getFullYear()} ${meta.author.name}`;

export default [
  {
    input: "index.js",
    output: {
      extend: true,
      //banner: copyright,
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
                //preamble: copyright
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
                //preamble: copyright
            }
        })
    ],
    output: {
      extend: true,
      file: "dist/druid.esm.js",
      format: "es",
      indent: false,
      name: "druid"
    }
  }
];