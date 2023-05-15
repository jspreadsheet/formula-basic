const { expect } = require('chai');

describe('Perform the function', () => {
    it('TODAY with simple arguments successfully', () => {
        const today = new Date()
        today.setMinutes(today.getMinutes() - today.getTimezoneOffset())

        expect(formula(`TODAY()`).getDate()).to.eq(today.getDate())
    })
});
