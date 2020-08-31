import resolve from "@rollup/plugin-node-resolve";
import { terser } from "rollup-plugin-terser";
import json from '@rollup/plugin-json';
import * as meta from "./package.json";
import jsdoc from 'rollup-plugin-jsdoc';

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
        resolve({
            customResolveOptions: {
                moduleDirectory: 'node_modules'
              }
        }),
        jsdoc({
            args: ["-r", "-d", "docs"],
            config: "jsdoc.config.json",
        }),
        json()
    ]
  },
  {
    input: "index.js",
    plugins: [
        resolve({
            customResolveOptions: {
                moduleDirectory: 'node_modules'
              }
        }),
        terser({
            format: {
                preamble: copyright
            }}),
        json()
    ],
    output: {
      extend: true,
      file: "dist/druid.min.js",
      format: "umd",
      indent: false,
      name: "druid"
    }
  }
];