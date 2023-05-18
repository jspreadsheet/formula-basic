const path = require('path');

module.exports = {
    target: ['web', 'es5'],
    entry: './src/formula.js',
    mode: 'production',
    output: {
        filename: 'formula.js',
        library: 'Formula'
    },
    optimization: {
        minimize: true
    },
    devServer: {
        static : {
            directory : path.join(__dirname, "/dist")
        },
        headers: {
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Methods": "GET, POST, PUT, DELETE, PATCH, OPTIONS",
            "Access-Control-Allow-Headers": "X-Requested-With, content-type, Authorization"
        },
        port: 3007,
        devMiddleware: {
            publicPath: "https://localhost:3000/dist/",
        },
        hot: "only",
    },
    stats: { warnings:false },
};
