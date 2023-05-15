#! /usr/bin/env node

let formula = require("./src/formula.js")

global.formula = formula;

exports.mochaHooks = {
    afterEach(done) {
        // destroy
        done();
    },
};
