const { expect } = require('chai')

describe('Perform the function', () => {
    it('AVERAGE with simple arguments successfully', () => {
        expect(formula('AVERAGE()')).to.be.throw
        expect(formula('AVERAGE(1)')).to.eq(1)
        expect(formula('AVERAGE(1, 1)')).to.eq(1)
        expect(formula('AVERAGE(-1, 1)')).to.eq(0)
        expect(formula('AVERAGE(-1, -1)')).to.eq(-1)
        expect(formula('AVERAGE(5, 5, 5, 5)')).to.eq(5)
        expect(formula('AVERAGE(1, 3, 5, 7)')).to.eq(4)
        expect(formula('AVERAGE(-10000, 10000, 345, -345)')).to.eq(0)
        expect(formula('AVERAGE(2.5, 8, 22, 25.7)')).to.be.approximately(14.55, 0.0001)
    })

    it('AVERAGE with named arguments successfully', () => {
        expect(formula('AVERAGE(A1, A2)', { A1: 4, A2: 9 })).to.eq(6.5)
        expect(formula('AVERAGE(A1, A2, B8)', { A1: 4, A2: 9, B8: -22 })).to.eq(-3)
        expect(formula('AVERAGE(A1, A2, B8)', { A1: 4, A2: 9, B8: -14 })).to.be.approximately(-0.3333, 0.0001)
    })

    it('MAX with simple arguments successfully', () => {
        expect(formula('MAX()')).to.be.throw
        expect(formula('MAX(1, 5)')).to.eq(5)
        expect(formula('MAX(7, 1, 2, 3, 4, 5, 6)')).to.eq(7)
        expect(formula('MAX(-1, 1)')).to.eq(1)
        expect(formula('MAX(-100, -1)')).to.eq(-1)
        expect(formula('MAX(1, 3, 5, 7)')).to.eq(7)
        expect(formula('MAX(-10000, 10000, 345, -345)')).to.eq(10000)
        expect(formula('MAX(2.5, 8, 22, 25.7)')).to.eq(25.7)
    })

    it('MAX with named arguments successfully', () => {
        expect(formula('MAX(A1, A2)', { A1: 4, A2: 9 })).to.eq(9)
        expect(formula('MAX(A1, A2, B8)', { A1: 4, A2: 9, B8: -22 })).to.eq(9)
        expect(formula('MAX(A1, A2, B8)', { A1: 4, A2: 9, B8: 45 })).to.eq(45)
        expect(formula('MAX(A1:A5)', { A1: 4, A2: 9, A3: 45, A4: 43, A5: -10 })).to.eq(45)
        expect(formula('MAX(A1:B2)', { A1: 4, A2: 9, B1: 3, B2: 11, B3: 1000 })).to.eq(11)
    })

    it('MIN with simple arguments successfully', () => {
        expect(formula('MIN()')).to.be.throw
        expect(formula('MIN(1, 5)')).to.eq(1)
        expect(formula('MIN(7, 1, 2, 3, 4, 5, 6)')).to.eq(1)
        expect(formula('MIN(-1, 1)')).to.eq(-1)
        expect(formula('MIN(-100, -1)')).to.eq(-100)
        expect(formula('MIN(1, 3, 5, 7)')).to.eq(1)
        expect(formula('MIN(-10000, 10000, 345, -345)')).to.eq(-10000)
        expect(formula('MIN(2.5, 8, 22, 25.7)')).to.eq(2.5)
    })

    it('MIN with named arguments successfully', () => {
        expect(formula('MIN(A1, A2)', { A1: 4, A2: 9 })).to.eq(4)
        expect(formula('MIN(A1, A2, B8)', { A1: 4, A2: 9, B8: -22 })).to.eq(-22)
        expect(formula('MIN(A1, A2, B8)', { A1: 4, A2: 9, B8: 45 })).to.eq(4)
        expect(formula('MIN(A1:A5)', { A1: 4, A2: 9, A3: 45, A4: 43, A5: -10 })).to.eq(-10)
        expect(formula('MIN(A1:B2)', { A1: 4, A2: 9, B1: 3, B2: 11, B3: -1000 })).to.eq(3)
    })

    it('COUNT with simple arguments successfully', () => {
        expect(formula('COUNT()')).to.eq(0)
        expect(formula('COUNT(1, 5)')).to.eq(2)
        expect(formula('COUNT(7, 1, 2, 3, 4, 5, 6)')).to.eq(7)
    })

    it('COUNT with named arguments successfully', () => {
        expect(formula('COUNT(A1, A2)', { A1: 4, A2: 9 })).to.eq(2)
        expect(formula('COUNT(A1, A2, B8)', { A1: 4, A2: 9, B8: 45 })).to.eq(3)
        expect(formula('COUNT(A1:A5)', { A1: 4, A2: 9, A3: 45, A4: 43, A5: -10 })).to.eq(5)
        expect(formula('COUNT(A1:B2)', { A1: 4, A2: 9, B1: 3, B2: 11, B3: -1000 })).to.eq(4)
        expect(formula('COUNT(A1:B2)', { A1: 4, A2: 9, B1: null, B2: 12 })).to.eq(3)
    })

    it('COUNTIF with simple arguments successfully', () => {
        expect(formula(`COUNTIF(5, '>3')`)).to.eq(1)
        expect(formula(`COUNTIF(2, '>3')`)).to.eq(0)
    })

    it('COUNTIF with named arguments successfully', () => {
        expect(formula(`COUNTIF(A1:A5, '>5')`, { A1: 4, A2: 9, A3: 5, A4: 6, A5: 22 })).to.eq(3)
        expect(formula(`COUNTIF(A1:A5, '>4')`, { A1: 4, A2: 9, A3: 5, A4: 6, A5: 22 })).to.eq(4)
        expect(formula(`COUNTIF(A1:A5, '>3')`, { A1: 4, A2: 9, A3: 5, A4: 6, A5: 22 })).to.eq(5)
        expect(formula(`COUNTIF(A1:A5, '<3')`, { A1: 4, A2: 9, A3: 5, A4: 6, A5: 22 })).to.eq(0)
        expect(formula(`COUNTIF(A1:A5, '<22')`, { A1: 4, A2: 9, A3: 5, A4: 6, A5: 22 })).to.eq(4)
    })
})
