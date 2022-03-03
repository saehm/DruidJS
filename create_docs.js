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
    "plugin": ["@saehrimnir/dmd-plugin"],
    "seperator": true,
    "package": "./package.json",
    "module-index-format": "table"
}

jsdoc2md.clear();

const json = jsdoc2md.getJsdocDataSync(Object.assign({
    "source": {
        "include": "**/*.js",
        "exclude": ["node_modules", "test", "docs", ".*"]
    }
}));

console.log(JSON.stringify(json, null, 4))
        

async function create_docs() {
    for (const module of modules) {
        const config = Object.assign({}, default_config)
        config.files = `${module}/index.js`;
        config.output = `${module}/README.md`;
        //jsdoc2md.getNamepaths(config)
        /* const jsdoc_data = await jsdoc2md.getJsdocData(config);
        const json = await jsdoc2md.getTemplateData({source: JSON.stringify(jsdoc_data)})
        json.forEach(entry => {
            if (entry?.meta?.code?.name) {
                console.log(entry)
                entry.meta.code.name =  entry.meta.code.name.replace(/exports\./, "druid.")
                console.log(entry)
            }
        })
        config.data = json */
        const output = await jsdoc2md.render(config)
        fs.writeFileSync(config.output, output);
    }
    
    await jsdoc2md.clear();
}


create_docs()