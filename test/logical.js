const { expect } = require('chai');

describe('Perform the function', () => {
    it('IF with simple arguments successfully', () => {
        expect(formula(`IF(2 > 1, 'YES', 'NO')`)).to.eq('YES')
        expect(formula(`IF(0 > 1, 'YES', 'NO')`)).to.eq('NO')
        expect(formula(`IF(1 = 1, 'YES', 'NO')`)).to.eq('YES')
        expect(formula(`IF(0 = 1, 'YES', 'NO')`)).to.eq('NO')
    })

    xit('IF with named arguments successfully', () => {
        expect(formula(`IF(A1 > 1, 'YES', 'NO')`, { A1: 2 })).to.eq('YES')
        expect(formula(`IF(A1 > 1, 'YES', 'NO')`, { A1: 0 })).to.eq('NO')
        expect(formula(`IF(A1 > A2, 'YES', 'NO')`, { A1: 7567567, A2: 1 })).to.eq('YES')
        expect(formula(`IF(A1 > A2, 'YES', 'NO')`, { A1: 0, A2: 23521 })).to.eq('NO')
        expect(formula(`IF(A1 > A2, D1, D2)`, { A1: 7567567, A2: 1, D1: 22, D2: -22 })).to.eq(22)
        expect(formula(`IF(A1 > A2, D1, D2)`, { A1: 7567567, A2: 1, D1: 'VERDADEIRO', D2: 'FALSO' })).to.eq('VERDADEIRO')
        expect(formula(`IF(A1 > A2, D1, D2)`, { A1: 0, A2: 23521, D1: 'VERDADEIRO', D2: 'FALSO' })).to.eq('FALSO')
    })
});
