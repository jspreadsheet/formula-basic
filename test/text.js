const { expect } = require('chai')

describe('Perform the function', () => {
    it('AVERAGE with simple arguments successfully', () => {
        expect(formula(`CONCATENATE('Hello', 'World')`)).to.eq('HELLOWORLD')
        expect(formula(`CONCATENATE('The', ' Lost', ' Symbol')`)).to.eq('THE LOST SYMBOL')
        expect(formula(`CONCATENATE('The', ' Lost', ' Symbol')`)).to.eq('THE LOST SYMBOL')
        expect(formula(`CONCATENATE(1, 2, 3)`)).to.eq('123')
    })

    xit('AVERAGE with named arguments successfully', () => {
        expect(formula('CONCATENATE(A1, A2)', { A1: 'HELLO', A2: 'WORLD' })).to.eq('HELLOWORLD')
    })
})
