# Formula

Formula Basic is an expression evaluator designed to handle Excel Formulas. It utilizes the [formulajs/formulajs](https://github.com/formulajs/formulajs) library function to process and calculate these Formulas.

## Installation

### npm

To install in your project using npm, run the following command:

```bash
$ npm i @jspreadsheet/formula-basic
```

### CDN

To use Formula via a CDN, include the following script tags in your HTML file:

```html
<script src="https://jsuites.net/v4/jsuites.js"></script>
<script src="https://cdn.jsdelivr.net/npm/@formulajs/formulajs/lib/browser/formula.min.js"></script>
<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/@jspreadsheet/formula/dist/index.min.js"></script>
```

## Usage

After installation, Formula can be utilized in your code by requiring it into your JavaScript file.

```javascript
const formula = require('@jspreadsheet/formula')
```

Doing basic math operations:

```javascript
formula('1 + 1')
// result -> 2
formula('2 * (4 + 4)')
// result -> 16
formula('10/(1+1)')
// result -> 5
formula('10**2')
// result -> 100
```

Using excel formulas:

```javascript
formula('SUM(2, 2, 4)')
// result -> 8
formula('AVERAGE(50, 55, 60)')
// result -> 55
formula('TODAY()')
```

Variables can be defined to use in the second argument of the function, for example:

```javascript
// Random variable
formula('HEYYOU', { HEYYOU: 5 })
// result -> 5
formula('HELLO + WORLD', { HELLO: 2, WORLD: -4 })
// result -> -2

// Excel cell names
formula('A1 + A2', { A1: 2, A2: 2 })
// result -> 4
formula('SUM(A1, A2, A3)', { A1: 6, A2: 5, A3: 10 })
// result -> 21
formula('COUNT(A1:B2)', { A1: 1, A2: 1, B1: 1, B2: null })
// result -> 3
```

The syntax **A1:A5** generates an array consisting of the values found by the formula for each A in the range from 1 to 5, in this case: [A1, A2, A3, A4, A5].

## Development

### Running the project

To run the project in development mode, use the following commands:

```bash
$ npm i
$ npm start
```

This will start a web-server with formula available in console, as a playground.

### Running Tests

After installing the packages run:

```bash
$ npm run test
```

To see more details in a browser:

```bash
$ npm run test:browser
```

To have more information about test coverage:

```bash
$ npm run test:coverage
```

## Contributing

Formula basic is an open source project and contributions are welcome! If you find a bug or have a feature request, please open an issue on GitHub. If you'd like to contribute code, please fork the repository and submit a pull request.

Ensure that you run the formatting plugins to maintain consistent code patterns. You can use the following command to do that:

```bash
$ npm run format
```


## License

Formula basic is released under the MIT.
