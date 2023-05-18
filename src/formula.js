/**
 * Jspreadsheet Extensions
 * Extension: Formula Basic
 * License: This is a free software MIT
 *
 * https://jspreadsheet.com
 */

import * as formulajs from '@formulajs/formulajs'

const Formula = function (scope) {
    // Based on sutoiku work (https://github.com/sutoiku)

    function getValue(obj, path) {
        const keys = path.split('.')
        let current = obj

        for (const key of keys) {
            if (current === undefined || current === null) {
                return undefined
            }
            current = current[key]
        }

        return current
    }

    for (let i = 0; i < Object.keys(formulajs).length; i++) {
        let method = Object.keys(formulajs)[i]
        let keys = []
        let values
        if (typeof formulajs[method] == 'object') {
            keys = Object.keys(formulajs[method])
            values = Object.values(formulajs[method])
            for (let a = 0; a < values.length; a++) {
                if (typeof values[a] == 'object') {
                    let subMethod = keys[a]
                    if (formulajs[method][subMethod]) {
                        keys = [
                            ...keys,
                            ...Object.keys(formulajs[method][subMethod]).map((a) => subMethod + '.' + a)
                        ] // Line too heavy, need refactor
                        keys.splice(keys.indexOf(subMethod), 1)
                    }
                }
            }
        }

        if (keys.length < 1) {
            scope[method] = formulajs[method]
        } else {
            for (let j = 0; j < keys.length; j++) {
                if (typeof getValue(formulajs[method], keys[j]) == 'function') {
                    scope[method] = getValue(formulajs[method], keys[j])
                }
            }
        }
    }

    /**
     * Instance execution helpers
     */
    var x = null
    var y = null
    var instance = null

    scope['TABLE'] = function () {
        return instance
    }
    scope['COLUMN'] = scope['COL'] = function () {
        if (instance.tracking) {
            instance.tracking.push(F.getColumnNameFromCoords(parseInt(x), parseInt(y)))
        }

        return parseInt(x) + 1
    }
    scope['ROW'] = function () {
        if (instance.tracking) {
            instance.tracking.push(F.getColumnNameFromCoords(parseInt(x), parseInt(y)))
        }

        return parseInt(y) + 1
    }
    scope['CELL'] = function () {
        return F.getColumnNameFromCoords(x, y)
    }
    scope['VALUE'] = function (col, row, processed) {
        return instance.getValueFromCoords(parseInt(col) - 1, parseInt(row) - 1, processed)
    }
    scope['THISROWCELL'] = function (col) {
        return instance.getValueFromCoords(parseInt(col) - 1, parseInt(y))
    }

    // Secure formula
    var secureFormula = function (oldValue, runtime) {
        var newValue = ''
        var inside = 0

        var special = ['=', '!', '>', '<']

        for (var i = 0; i < oldValue.length; i++) {
            if (oldValue[i] == '"') {
                if (inside == 0) {
                    inside = 1
                } else {
                    inside = 0
                }
            }

            if (inside == 1) {
                newValue += oldValue[i]
            } else {
                newValue += oldValue[i].toUpperCase()

                if (runtime == true) {
                    if (
                        i > 0 &&
                        oldValue[i] == '=' &&
                        special.indexOf(oldValue[i - 1]) == -1 &&
                        special.indexOf(oldValue[i + 1]) == -1
                    ) {
                        newValue += '='
                    }
                }
            }
        }

        // Adapt to JS
        newValue = newValue.replace(/\^/g, '**')
        newValue = newValue.replace(/\<\>/g, '!=')
        newValue = newValue.replace(/\&/g, '+')
        newValue = newValue.replace(/\$/g, '')

        return newValue
    }

    // Convert range tokens
    var tokensUpdate = function (tokens, e) {
        for (var index = 0; index < tokens.length; index++) {
            var f = F.getTokensFromRange(tokens[index])
            e = e.replace(tokens[index], '[' + f.join(',') + ']')
        }
        return e
    }

    var isNumeric = function (num) {
        if (typeof (num) === 'string') {
            num = num.trim();
        }
        return !isNaN(num) && num !== null && num !== '';
    }

    var F = function (expression, variables, i, j, obj) {
        // Global helpers
        instance = obj
        x = i
        y = j
        // String
        var s = ''
        var parent = {}
        if (variables) {
            if (variables.size) {
                tokens = null
                variables.forEach(function (v, k) {
                    // Replace ! per dot
                    t = k.replace(/!/g, '.')
                    // Exists parent
                    if (t.indexOf('.') !== -1) {
                        t = t.split('.')
                        parent[t[0]] = true
                    }
                })
                t = Object.keys(parent)
                for (let i = 0; i < t.length; i++) {
                    s += 'var ' + t[i] + ' = {};'
                }

                variables.forEach(function (v, k) {
                    // Replace ! per dot
                    t = k.replace(/!/g, '.')
                    if (v !== null && !isNumeric(v)) {
                        tokens = v.match(/(('.*?'!)|(\w*!))?(\$?[A-Z]+\$?[0-9]*):(\$?[A-Z]+\$?[0-9]*)?/g)
                        if (tokens && tokens.length) {
                            v = updateRanges(tokens, v)
                        }
                    }

                    if (t.indexOf('.') > 0) {
                        s += t + ' = ' + variables.get(k) + ';\n'
                    } else {
                        s += 'var ' + t + ' = ' + v + ';\n'
                    }
                })
            } else {
                var keys = Object.keys(variables)
                if (keys.length) {
                    var parent = {}
                    for (var i = 0; i < keys.length; i++) {
                        // Replace ! per dot
                        t = keys[i].replace(/\!/g, '.')
                        // Exists parent
                        if (t.indexOf('.') > 0) {
                            var t = t.split('.')
                            parent[t[0]] = {}
                        }
                    }
                    var t = Object.keys(parent)
                    for (var i = 0; i < t.length; i++) {
                        s += 'var ' + t[i] + ' = {};'
                    }

                    for (var i = 0; i < keys.length; i++) {
                        // Replace ! per dot
                        var t = keys[i].replace(/\!/g, '.')

                        // Update range
                        if (variables[keys[i]] !== null && !isNumeric(variables[keys[i]])) {
                            var tokens = variables[keys[i]].match(
                                /(('.*?'!)|(\w*!))?(\$?[A-Z]+\$?[0-9]*):(\$?[A-Z]+\$?[0-9]*)?/g
                            )
                            if (tokens && tokens.length) {
                                variables[keys[i]] = tokensUpdate(tokens, variables[keys[i]])
                            }
                        }

                        if (t.indexOf('.') > 0) {
                            s += t + " = " + variables[keys[i]] + ";\n"
                        } else {
                            s += "var " + t + " = " + variables[keys[i]] + ";\n"
                        }
                    }
                }
            }
        }
        // Remove $
        expression = expression.replace(/\$/g, '')
        // Replace ! per dot
        expression = expression.replace(/\!/g, '.')
        // Adapt to JS
        expression = secureFormula(expression, true)
        // Update range
        var tokens = expression.match(/(('.*?'!)|(\w*!))?(\$?[A-Z]+\$?[0-9]*):(\$?[A-Z]+\$?[0-9]*)?/g)
        if (tokens && tokens.length) {
            expression = tokensUpdate(tokens, expression)
        }
        // Calculate
        var result = new Function(s + "; return " + expression)()
        if (result === null) {
            result = 0
        }

        return result
    }

    F.setFormula = function (o) {
        var k = Object.keys(o)
        for (var i = 0; i < k.length; i++) {
            if (typeof o[k[i]] == 'function') {
                scope[k[i]] = o[k[i]]
            }
        }
    }

    F.basic = true

    return F
}

let formula = Formula(typeof window === 'undefined' ? global : window)
export default formula;