//import ascii from "rollup-plugin-ascii";
import resolve from "@rollup/plugin-node-resolve";
import {terser} from "rollup-plugin-terser";
//import {eslint} from 'rollup-plugin-eslint';
import * as meta from "./package.json";
import jsdoc from 'rollup-plugin-jsdoc';
//import babel from "rollup-plugin-babel";

const copyright = `// ${meta.homepage} v${meta.version} Copyright ${(new Date).getFullYear()} ${meta.author.name}`;

export default [
  {
    input: "index",
    output: {
      extend: true,
      banner: copyright,
      file: "dist/druid.js",
      format: "umd",
      indent: false,
      name: "DruidJS"
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
      //eslint()
    ]
  },
  {
    input: "index",
    plugins: [
      terser({output: {preamble: copyright}}),
      //eslint()
    ],
    output: {
      extend: true,
      file: "dist/druid.min.js",
      format: "umd",
      indent: false,
      name: "DruidJS"
    }
  }
];