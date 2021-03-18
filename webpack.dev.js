// eslint-disable-next-line unicorn/prevent-abbreviations
const path = require('path');

module.exports = {
    entry: path.resolve(__dirname, './index.js'),
    output: {
        path: path.resolve(__dirname, 'dist'),
        filename: 'druid.js',
        library: 'druid',
        libraryTarget: 'umd',
        globalObject: 'this',
    },
    module: {
        rules: [
            {
                test: /\.(js)$/,
                exclude: /node_modules/,
                enforce: 'pre',
                use: ['source-map-loader', 'babel-loader'],
            },
        ],
    },
    devtool: 'eval-source-map',
    mode: 'development',
};
