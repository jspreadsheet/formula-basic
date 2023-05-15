const { expect } = require('chai');

describe('Perform a basic operation', () => {
    it('successfully', () => {
        expect(formula('1 + 1')).to.eq(2)
        expect(formula('2 + 5')).to.eq(7)
        expect(formula('8986 + 9823')).to.eq(18809)

        expect(formula('1 - 1')).to.eq(0)
        expect(formula('(-1) + (-1)')).to.eq(-2)
        expect(formula('0 - 999')).to.eq(-999)

        expect(formula('10 / 2')).to.eq(5)
        expect(formula('60 / 1')).to.eq(60)
        expect(formula('-10 / 2')).to.eq(-5)
        expect(formula('5 / 0')).to.eq(Infinity)
        expect(formula('-5 / 0')).to.eq(-Infinity)

        expect(formula('5 * 5')).to.eq(25)
        expect(formula('5 * 5 + 1')).to.eq(26)
        expect(formula('5 * (5 + 1)')).to.eq(30)
        expect(formula('5 * 5 + 10 / 2')).to.eq(30)

        expect(formula('2 ** 3')).to.eq(8)
        expect(formula('10 ** 0')).to.eq(1)
        expect(formula('(-2) ** 2')).to.eq(4)

        expect(formula('10 % 3')).to.eq(1)
        expect(formula('5 % 2')).to.eq(1)
        expect(formula('(-10) % 3')).to.eq(-1)
    })
});
