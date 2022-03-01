import jsdoc2md from "jsdoc-to-markdown"
import fs from "fs";

const modules = [
    "clustering",
    "datastructure",
    "dimred",
    "knn",
    "linear_algebra",
    "matrix",
    "metrics",
    "numerical",
    "optimization",
    "util"
]

const default_config = {
    "files": "*.js",
    "param-list-format": "list",
    "plugin": ["dmd-clear"],
    "seperator": true,
    
}

for (const module of modules) {
    const config = Object.assign({}, default_config)
    config.files = `${module}/*.js`;
    config.output = `${module}/README.md`;
    jsdoc2md.render(config).then(output => {
        fs.writeFileSync(config.output, output);
    })
}
