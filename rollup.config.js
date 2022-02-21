import resolve from "@rollup/plugin-node-resolve";
import { terser } from "rollup-plugin-terser";
import jsdoc from 'rollup-plugin-jsdoc';
import json from '@rollup/plugin-json';
import meta from "./package.json";

const copyright = `// ${meta.homepage} v${meta.version} Copyright ${(new Date).getFullYear()} ${meta.author.name}`;

const onwarn = (message, warn) => {
  if (message.code === "CIRCULAR_DEPENDENCY") return;
  warn(message);
}

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
        json({
            compact: true,
            exclude: 'node_modules/**',
        }),
        resolve(),
        jsdoc({
            args: ["-r", "-d", "docs"],
            config: "jsdoc.config.json",
        }), 
    ],
    onwarn
  }, {
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
      sourcemap: 'inline',
      extend: true,
      file: "dist/druid.min.js",
      format: "umd",
      indent: false,
      name: "druid"
    },
    onwarn
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
      sourcemap: 'inline',
      extend: true,
      file: "dist/druid.esm.js",
      format: "es",
      indent: false,
      name: "druid"
    },
    onwarn
  }
];