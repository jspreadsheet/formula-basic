const path = require('path');

class MyPlugin {
    apply(compiler) {
        compiler.hooks.emit.tap('MyPlugin', (compilation) => {
            // Get the bundled file name
            const fileName = Object.keys(compilation.assets)[0];

            // Get the bundled file content
            const fileContent = compilation.assets[fileName].source();

            const header = `;(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
    typeof define === 'function' && define.amd ? define(factory) :
    global.formula = factory();
}(this, (function () {`;

            const footer = `    return Formula;
})));`;

            // Updated file content with custom content added
            const updatedFileContent = header + '\n\n' + fileContent + '\n\n' + footer;

            // Replace the bundled file content with updated content
            compilation.assets[fileName] = {
                source: () => updatedFileContent,
                size: () => updatedFileContent.length,
            };
        });
    }
}

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
    plugins: [
        new MyPlugin(),
    ],
    stats: { warnings:false },
};
