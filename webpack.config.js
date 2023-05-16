const path = require('path')

module.exports = (env, argv) => {
    const prod = {
        target: ['web', 'es5'],
        entry: {
            index: './src/index.js'
        },
        mode: 'production',
        output: {
            library: {
                name: 'formula-basic',
                type: 'umd',
                export: ['default']
            },
            globalObject: 'this',
            filename: '[name].js'
        },
        plugins: [],
        optimization: {
            minimize: false
        },
        stats: { warnings: false }
    }

    const dev = {
        target: 'web',
        entry: './src/index.js',
        mode: 'development',
        optimization: {
            minimize: false
        },
        output: {
            filename: 'index.js',
            path: path.resolve(__dirname, 'dist')
        },
        module: {
            rules: [
                {
                    test: /\.css$/,
                    use: ['style-loader', 'css-loader']
                }
            ]
        },
        devServer: {
            // contentBase
            static: {
                directory: path.join(__dirname, '/')
            },
            headers: {
                'Access-Control-Allow-Origin': '*',
                'Access-Control-Allow-Methods': 'GET, POST, PUT, DELETE, PATCH, OPTIONS',
                'Access-Control-Allow-Headers': 'X-Requested-With, content-type, Authorization'
            },
            port: 3005,
            devMiddleware: {
                publicPath: 'https://localhost:3000/dist/'
            },
            hot: 'only'
        },
        stats: { warnings: false }
    }

    return argv.mode === 'production' ? prod : dev
}
