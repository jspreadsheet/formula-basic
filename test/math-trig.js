const { expect } = require('chai');

describe('Perform the function', () => {
    it('SUM with simple arguments successfully', () => {
        expect(formula('SUM()')).to.eq(0)
        expect(formula('SUM(1)')).to.eq(1)
        expect(formula('SUM(1, 1)')).to.eq(2)
        expect(formula('SUM(-1, 1)')).to.eq(0)
        expect(formula('SUM(-1, -1)')).to.eq(-2)
        expect(formula('SUM(5, 5, 5, 5)')).to.eq(20)
    })

    it('SUM with named arguments successfully', () => {
        expect(formula('SUM(A1, A2)', { A1: 4, A2: 9 })).to.eq(13)
        expect(formula('SUM(A1, A2, B8)', { A1: 4, A2: 9, B8: 7 })).to.eq(20)
        expect(formula('SUM(A1, A2, B8)', { A1: 4, A2: 9, B8: -14 })).to.eq(-1)
    })

    it('SUMIF with named arguments successfully', () => {
        expect(formula(`SUMIF(A1:A2, '>3')`, { A1: 4, A2: 9 })).to.eq(13)
        expect(formula(`SUMIF(A1:A5, '>5')`, { A1: 4, A2: 9, A3: 3, A4: 2, A5: 7 })).to.eq(16)
        expect(formula(`SUMIF(A1:A5, '<5')`, { A1: 4, A2: 9, A3: 3, A4: 2, A5: 7 })).to.eq(9)
        expect(formula(`SUMIF(A1:B2, '<0')`, { A1: -4, A2: -5, B1: -12, B2: -3, A5: -100 })).to.eq(-24)
    })

    it('SUMPRODUCT with named arguments successfully', () => {
        expect(formula(`SUMPRODUCT(A1:A2, B1:B2)`, { A1: 4, A2: 9, B1: 4, B2: 4 })).to.eq(52)
        expect(formula(`SUMPRODUCT(A1:A3, B1:B3)`, { A1: 4, A2: 9, A3: 1, B1: 4, B2: 4, B3: -10 })).to.eq(42)
        expect(formula(`SUMPRODUCT(A1:A5, B1:B5)`, { A1: 4, A2: 9, A3: 1, A4: 15, A5: -3, B1: 4, B2: 4, B3: -10, B4: 0, B5: -3 })).to.eq(51)
    })
});
