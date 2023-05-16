(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory();
	else if(typeof define === 'function' && define.amd)
		define([], factory);
	else if(typeof exports === 'object')
		exports["formula-basic"] = factory();
	else
		root["formula-basic"] = factory();
})(this, function() {
return /******/ (function() { // webpackBootstrap
/******/ 	var __webpack_modules__ = ({

/***/ 765:
/***/ (function(__unused_webpack_module, exports) {

/* bessel.js (C) 2013-present SheetJS -- http://sheetjs.com */
/* vim: set ts=2: */
/*exported BESSEL */
var BESSEL;
(function (factory) {
  /*jshint ignore:start */
  if(typeof DO_NOT_EXPORT_BESSEL === 'undefined') {
    if(true) {
      factory(exports);
    } else {}
  } else {
    factory(BESSEL = {});
  }
  /*jshint ignore:end */
}(function(BESSEL) {
BESSEL.version = '1.0.2';
var M = Math;

function _horner(arr, v) { for(var i = 0, z = 0; i < arr.length; ++i) z = v * z + arr[i]; return z; }
function _bessel_iter(x, n, f0, f1, sign) {
  if(n === 0) return f0;
  if(n === 1) return f1;
  var tdx = 2 / x, f2 = f1;
  for(var o = 1; o < n; ++o) {
    f2 = f1 * o * tdx + sign * f0;
    f0 = f1; f1 = f2;
  }
  return f2;
}
function _bessel_wrap(bessel0, bessel1, name, nonzero, sign) {
  return function bessel(x,n) {
    if(nonzero) {
      if(x === 0) return (nonzero == 1 ? -Infinity : Infinity);
      else if(x < 0) return NaN;
    }
    if(n === 0) return bessel0(x);
    if(n === 1) return bessel1(x);
    if(n < 0) return NaN;
    n|=0;
    var b0 = bessel0(x), b1 = bessel1(x);
    return _bessel_iter(x, n, b0, b1, sign);
  };
}
var besselj = (function() {
  var W = 0.636619772; // 2 / Math.PI

  var b0_a1a = [57568490574.0, -13362590354.0, 651619640.7, -11214424.18, 77392.33017, -184.9052456].reverse();
  var b0_a2a = [57568490411.0, 1029532985.0, 9494680.718, 59272.64853, 267.8532712, 1.0].reverse();
  var b0_a1b = [1.0, -0.1098628627e-2, 0.2734510407e-4, -0.2073370639e-5, 0.2093887211e-6].reverse();
  var b0_a2b = [-0.1562499995e-1, 0.1430488765e-3, -0.6911147651e-5, 0.7621095161e-6, -0.934935152e-7].reverse();

  function bessel0(x) {
    var a=0, a1=0, a2=0, y = x * x;
    if(x < 8) {
      a1 = _horner(b0_a1a, y);
      a2 = _horner(b0_a2a, y);
      a = a1 / a2;
    } else {
      var xx = x - 0.785398164;
      y = 64 / y;
      a1 = _horner(b0_a1b, y);
      a2 = _horner(b0_a2b, y);
      a = M.sqrt(W/x)*(M.cos(xx)*a1-M.sin(xx)*a2*8/x);
    }
    return a;
  }

  var b1_a1a = [72362614232.0, -7895059235.0, 242396853.1, -2972611.439, 15704.48260, -30.16036606].reverse();
  var b1_a2a = [144725228442.0, 2300535178.0, 18583304.74, 99447.43394, 376.9991397, 1.0].reverse();
  var b1_a1b = [1.0, 0.183105e-2, -0.3516396496e-4, 0.2457520174e-5, -0.240337019e-6].reverse();
  var b1_a2b = [0.04687499995, -0.2002690873e-3, 0.8449199096e-5, -0.88228987e-6, 0.105787412e-6].reverse();

  function bessel1(x) {
    var a=0, a1=0, a2=0, y = x*x, xx = M.abs(x) - 2.356194491;
    if(Math.abs(x)< 8) {
      a1 = x*_horner(b1_a1a, y);
      a2 = _horner(b1_a2a, y);
      a = a1 / a2;
    } else {
      y = 64 / y;
      a1=_horner(b1_a1b, y);
      a2=_horner(b1_a2b, y);
      a=M.sqrt(W/M.abs(x))*(M.cos(xx)*a1-M.sin(xx)*a2*8/M.abs(x));
      if(x < 0) a = -a;
    }
    return a;
  }

  return function besselj(x, n) {
    n = Math.round(n);
    if(!isFinite(x)) return isNaN(x) ? x : 0;
    if(n < 0) return ((n%2)?-1:1)*besselj(x, -n);
    if(x < 0) return ((n%2)?-1:1)*besselj(-x, n);
    if(n === 0) return bessel0(x);
    if(n === 1) return bessel1(x);
    if(x === 0) return 0;

    var ret=0.0;
    if(x > n) {
      ret = _bessel_iter(x, n, bessel0(x), bessel1(x),-1);
    } else {
      var m=2*M.floor((n+M.floor(M.sqrt(40*n)))/2);
      var jsum=false;
      var bjp=0.0, sum=0.0;
      var bj=1.0, bjm = 0.0;
      var tox = 2 / x;
      for (var j=m;j>0;j--) {
        bjm=j*tox*bj-bjp;
        bjp=bj;
        bj=bjm;
        if (M.abs(bj) > 1E10) {
          bj *= 1E-10;
          bjp *= 1E-10;
          ret *= 1E-10;
          sum *= 1E-10;
        }
        if (jsum) sum += bj;
        jsum=!jsum;
        if (j == n) ret=bjp;
      }
      sum=2.0*sum-bj;
      ret /= sum;
    }
    return ret;
  };
})();
var bessely = (function() {
  var W = 0.636619772;

  var b0_a1a = [-2957821389.0, 7062834065.0, -512359803.6, 10879881.29, -86327.92757, 228.4622733].reverse();
  var b0_a2a = [40076544269.0, 745249964.8, 7189466.438, 47447.26470, 226.1030244, 1.0].reverse();
  var b0_a1b = [1.0, -0.1098628627e-2, 0.2734510407e-4, -0.2073370639e-5, 0.2093887211e-6].reverse();
  var b0_a2b = [-0.1562499995e-1, 0.1430488765e-3, -0.6911147651e-5, 0.7621095161e-6, -0.934945152e-7].reverse();

  function bessel0(x) {
    var a=0, a1=0, a2=0, y = x * x, xx = x - 0.785398164;
    if(x < 8) {
      a1 = _horner(b0_a1a, y);
      a2 = _horner(b0_a2a, y);
      a = a1/a2 + W * besselj(x,0) * M.log(x);
    } else {
      y = 64 / y;
      a1 = _horner(b0_a1b, y);
      a2 = _horner(b0_a2b, y);
      a = M.sqrt(W/x)*(M.sin(xx)*a1+M.cos(xx)*a2*8/x);
    }
    return a;
  }

  var b1_a1a = [-0.4900604943e13, 0.1275274390e13, -0.5153438139e11, 0.7349264551e9, -0.4237922726e7, 0.8511937935e4].reverse();
  var b1_a2a = [0.2499580570e14, 0.4244419664e12, 0.3733650367e10, 0.2245904002e8, 0.1020426050e6, 0.3549632885e3, 1].reverse();
  var b1_a1b = [1.0, 0.183105e-2, -0.3516396496e-4, 0.2457520174e-5, -0.240337019e-6].reverse();
  var b1_a2b = [0.04687499995, -0.2002690873e-3, 0.8449199096e-5, -0.88228987e-6, 0.105787412e-6].reverse();

  function bessel1(x) {
    var a=0, a1=0, a2=0, y = x*x, xx = x - 2.356194491;
    if(x < 8) {
      a1 = x*_horner(b1_a1a, y);
      a2 = _horner(b1_a2a, y);
      a = a1/a2 + W * (besselj(x,1) * M.log(x) - 1 / x);
    } else {
      y = 64 / y;
      a1=_horner(b1_a1b, y);
      a2=_horner(b1_a2b, y);
      a=M.sqrt(W/x)*(M.sin(xx)*a1+M.cos(xx)*a2*8/x);
    }
    return a;
  }

  return _bessel_wrap(bessel0, bessel1, 'BESSELY', 1, -1);
})();
var besseli = (function() {
  var b0_a = [1.0, 3.5156229, 3.0899424, 1.2067492, 0.2659732, 0.360768e-1, 0.45813e-2].reverse();
  var b0_b = [0.39894228, 0.1328592e-1, 0.225319e-2, -0.157565e-2, 0.916281e-2, -0.2057706e-1, 0.2635537e-1, -0.1647633e-1, 0.392377e-2].reverse();

  function bessel0(x) {
    if(x <= 3.75) return _horner(b0_a, x*x/(3.75*3.75));
    return M.exp(M.abs(x))/M.sqrt(M.abs(x))*_horner(b0_b, 3.75/M.abs(x));
  }

  var b1_a = [0.5, 0.87890594, 0.51498869, 0.15084934, 0.2658733e-1, 0.301532e-2, 0.32411e-3].reverse();
  var b1_b = [0.39894228, -0.3988024e-1, -0.362018e-2, 0.163801e-2, -0.1031555e-1, 0.2282967e-1, -0.2895312e-1, 0.1787654e-1, -0.420059e-2].reverse();

  function bessel1(x) {
    if(x < 3.75) return x * _horner(b1_a, x*x/(3.75*3.75));
    return (x < 0 ? -1 : 1) * M.exp(M.abs(x))/M.sqrt(M.abs(x))*_horner(b1_b, 3.75/M.abs(x));
  }

  return function besseli(x, n) {
    n = Math.round(n);
    if(n === 0) return bessel0(x);
    if(n === 1) return bessel1(x);
    if(n < 0) return NaN;
    if(M.abs(x) === 0) return 0;
    if(x == Infinity) return Infinity;

    var ret = 0.0, j, tox = 2 / M.abs(x), bip = 0.0, bi=1.0, bim=0.0;
    var m=2*M.round((n+M.round(M.sqrt(40*n)))/2);
    for (j=m;j>0;j--) {
      bim=j*tox*bi + bip;
      bip=bi; bi=bim;
      if (M.abs(bi) > 1E10) {
        bi *= 1E-10;
        bip *= 1E-10;
        ret *= 1E-10;
      }
      if(j == n) ret = bip;
    }
    ret *= besseli(x, 0) / bi;
    return x < 0 && (n%2) ? -ret : ret;
  };

})();

var besselk = (function() {
  var b0_a = [-0.57721566, 0.42278420, 0.23069756, 0.3488590e-1, 0.262698e-2, 0.10750e-3, 0.74e-5].reverse();
  var b0_b = [1.25331414, -0.7832358e-1, 0.2189568e-1, -0.1062446e-1, 0.587872e-2, -0.251540e-2, 0.53208e-3].reverse();

  function bessel0(x) {
    if(x <= 2) return -M.log(x/2) * besseli(x,0) + _horner(b0_a, x*x/4);
    return M.exp(-x) / M.sqrt(x) * _horner(b0_b, 2/x);
  }

  var b1_a = [1.0, 0.15443144, -0.67278579, -0.18156897, -0.1919402e-1, -0.110404e-2, -0.4686e-4].reverse();
  var b1_b = [1.25331414, 0.23498619, -0.3655620e-1, 0.1504268e-1, -0.780353e-2, 0.325614e-2, -0.68245e-3].reverse();

  function bessel1(x) {
    if(x <= 2) return M.log(x/2) * besseli(x,1) + (1/x) * _horner(b1_a, x*x/4);
    return M.exp(-x)/M.sqrt(x)*_horner(b1_b, 2/x);
  }

  return _bessel_wrap(bessel0, bessel1, 'BESSELK', 2, 1);
})();
BESSEL.besselj = besselj;
BESSEL.bessely = bessely;
BESSEL.besseli = besseli;
BESSEL.besselk = besselk;
}));


/***/ }),

/***/ 162:
/***/ (function(module) {

(function (window, factory) {
    if (true) {
        module.exports = factory();
    } else {}
})(this, function () {
var jStat = (function(Math, undefined) {

// For quick reference.
var concat = Array.prototype.concat;
var slice = Array.prototype.slice;
var toString = Object.prototype.toString;

// Calculate correction for IEEE error
// TODO: This calculation can be improved.
function calcRdx(n, m) {
  var val = n > m ? n : m;
  return Math.pow(10,
                  17 - ~~(Math.log(((val > 0) ? val : -val)) * Math.LOG10E));
}


var isArray = Array.isArray || function isArray(arg) {
  return toString.call(arg) === '[object Array]';
};


function isFunction(arg) {
  return toString.call(arg) === '[object Function]';
}


function isNumber(num) {
  return (typeof num === 'number') ? num - num === 0 : false;
}


// Converts the jStat matrix to vector.
function toVector(arr) {
  return concat.apply([], arr);
}


// The one and only jStat constructor.
function jStat() {
  return new jStat._init(arguments);
}


// TODO: Remove after all references in src files have been removed.
jStat.fn = jStat.prototype;


// By separating the initializer from the constructor it's easier to handle
// always returning a new instance whether "new" was used or not.
jStat._init = function _init(args) {
  // If first argument is an array, must be vector or matrix.
  if (isArray(args[0])) {
    // Check if matrix.
    if (isArray(args[0][0])) {
      // See if a mapping function was also passed.
      if (isFunction(args[1]))
        args[0] = jStat.map(args[0], args[1]);
      // Iterate over each is faster than this.push.apply(this, args[0].
      for (var i = 0; i < args[0].length; i++)
        this[i] = args[0][i];
      this.length = args[0].length;

    // Otherwise must be a vector.
    } else {
      this[0] = isFunction(args[1]) ? jStat.map(args[0], args[1]) : args[0];
      this.length = 1;
    }

  // If first argument is number, assume creation of sequence.
  } else if (isNumber(args[0])) {
    this[0] = jStat.seq.apply(null, args);
    this.length = 1;

  // Handle case when jStat object is passed to jStat.
  } else if (args[0] instanceof jStat) {
    // Duplicate the object and pass it back.
    return jStat(args[0].toArray());

  // Unexpected argument value, return empty jStat object.
  // TODO: This is strange behavior. Shouldn't this throw or some such to let
  // the user know they had bad arguments?
  } else {
    this[0] = [];
    this.length = 1;
  }

  return this;
};
jStat._init.prototype = jStat.prototype;
jStat._init.constructor = jStat;


// Utility functions.
// TODO: for internal use only?
jStat.utils = {
  calcRdx: calcRdx,
  isArray: isArray,
  isFunction: isFunction,
  isNumber: isNumber,
  toVector: toVector
};


jStat._random_fn = Math.random;
jStat.setRandom = function setRandom(fn) {
  if (typeof fn !== 'function')
    throw new TypeError('fn is not a function');
  jStat._random_fn = fn;
};


// Easily extend the jStat object.
// TODO: is this seriously necessary?
jStat.extend = function extend(obj) {
  var i, j;

  if (arguments.length === 1) {
    for (j in obj)
      jStat[j] = obj[j];
    return this;
  }

  for (i = 1; i < arguments.length; i++) {
    for (j in arguments[i])
      obj[j] = arguments[i][j];
  }

  return obj;
};


// Returns the number of rows in the matrix.
jStat.rows = function rows(arr) {
  return arr.length || 1;
};


// Returns the number of columns in the matrix.
jStat.cols = function cols(arr) {
  return arr[0].length || 1;
};


// Returns the dimensions of the object { rows: i, cols: j }
jStat.dimensions = function dimensions(arr) {
  return {
    rows: jStat.rows(arr),
    cols: jStat.cols(arr)
  };
};


// Returns a specified row as a vector or return a sub matrix by pick some rows
jStat.row = function row(arr, index) {
  if (isArray(index)) {
    return index.map(function(i) {
      return jStat.row(arr, i);
    })
  }
  return arr[index];
};


// return row as array
// rowa([[1,2],[3,4]],0) -> [1,2]
jStat.rowa = function rowa(arr, i) {
  return jStat.row(arr, i);
};


// Returns the specified column as a vector or return a sub matrix by pick some
// columns
jStat.col = function col(arr, index) {
  if (isArray(index)) {
    var submat = jStat.arange(arr.length).map(function() {
      return new Array(index.length);
    });
    index.forEach(function(ind, i){
      jStat.arange(arr.length).forEach(function(j) {
        submat[j][i] = arr[j][ind];
      });
    });
    return submat;
  }
  var column = new Array(arr.length);
  for (var i = 0; i < arr.length; i++)
    column[i] = [arr[i][index]];
  return column;
};


// return column as array
// cola([[1,2],[3,4]],0) -> [1,3]
jStat.cola = function cola(arr, i) {
  return jStat.col(arr, i).map(function(a){ return a[0] });
};


// Returns the diagonal of the matrix
jStat.diag = function diag(arr) {
  var nrow = jStat.rows(arr);
  var res = new Array(nrow);
  for (var row = 0; row < nrow; row++)
    res[row] = [arr[row][row]];
  return res;
};


// Returns the anti-diagonal of the matrix
jStat.antidiag = function antidiag(arr) {
  var nrow = jStat.rows(arr) - 1;
  var res = new Array(nrow);
  for (var i = 0; nrow >= 0; nrow--, i++)
    res[i] = [arr[i][nrow]];
  return res;
};

// Transpose a matrix or array.
jStat.transpose = function transpose(arr) {
  var obj = [];
  var objArr, rows, cols, j, i;

  // Make sure arr is in matrix format.
  if (!isArray(arr[0]))
    arr = [arr];

  rows = arr.length;
  cols = arr[0].length;

  for (i = 0; i < cols; i++) {
    objArr = new Array(rows);
    for (j = 0; j < rows; j++)
      objArr[j] = arr[j][i];
    obj.push(objArr);
  }

  // If obj is vector, return only single array.
  return obj.length === 1 ? obj[0] : obj;
};


// Map a function to an array or array of arrays.
// "toAlter" is an internal variable.
jStat.map = function map(arr, func, toAlter) {
  var row, nrow, ncol, res, col;

  if (!isArray(arr[0]))
    arr = [arr];

  nrow = arr.length;
  ncol = arr[0].length;
  res = toAlter ? arr : new Array(nrow);

  for (row = 0; row < nrow; row++) {
    // if the row doesn't exist, create it
    if (!res[row])
      res[row] = new Array(ncol);
    for (col = 0; col < ncol; col++)
      res[row][col] = func(arr[row][col], row, col);
  }

  return res.length === 1 ? res[0] : res;
};


// Cumulatively combine the elements of an array or array of arrays using a function.
jStat.cumreduce = function cumreduce(arr, func, toAlter) {
  var row, nrow, ncol, res, col;

  if (!isArray(arr[0]))
    arr = [arr];

  nrow = arr.length;
  ncol = arr[0].length;
  res = toAlter ? arr : new Array(nrow);

  for (row = 0; row < nrow; row++) {
    // if the row doesn't exist, create it
    if (!res[row])
      res[row] = new Array(ncol);
    if (ncol > 0)
      res[row][0] = arr[row][0];
    for (col = 1; col < ncol; col++)
      res[row][col] = func(res[row][col-1], arr[row][col]);
  }
  return res.length === 1 ? res[0] : res;
};


// Destructively alter an array.
jStat.alter = function alter(arr, func) {
  return jStat.map(arr, func, true);
};


// Generate a rows x cols matrix according to the supplied function.
jStat.create = function  create(rows, cols, func) {
  var res = new Array(rows);
  var i, j;

  if (isFunction(cols)) {
    func = cols;
    cols = rows;
  }

  for (i = 0; i < rows; i++) {
    res[i] = new Array(cols);
    for (j = 0; j < cols; j++)
      res[i][j] = func(i, j);
  }

  return res;
};


function retZero() { return 0; }


// Generate a rows x cols matrix of zeros.
jStat.zeros = function zeros(rows, cols) {
  if (!isNumber(cols))
    cols = rows;
  return jStat.create(rows, cols, retZero);
};


function retOne() { return 1; }


// Generate a rows x cols matrix of ones.
jStat.ones = function ones(rows, cols) {
  if (!isNumber(cols))
    cols = rows;
  return jStat.create(rows, cols, retOne);
};


// Generate a rows x cols matrix of uniformly random numbers.
jStat.rand = function rand(rows, cols) {
  if (!isNumber(cols))
    cols = rows;
  return jStat.create(rows, cols, jStat._random_fn);
};


function retIdent(i, j) { return i === j ? 1 : 0; }


// Generate an identity matrix of size row x cols.
jStat.identity = function identity(rows, cols) {
  if (!isNumber(cols))
    cols = rows;
  return jStat.create(rows, cols, retIdent);
};


// Tests whether a matrix is symmetric
jStat.symmetric = function symmetric(arr) {
  var size = arr.length;
  var row, col;

  if (arr.length !== arr[0].length)
    return false;

  for (row = 0; row < size; row++) {
    for (col = 0; col < size; col++)
      if (arr[col][row] !== arr[row][col])
        return false;
  }

  return true;
};


// Set all values to zero.
jStat.clear = function clear(arr) {
  return jStat.alter(arr, retZero);
};


// Generate sequence.
jStat.seq = function seq(min, max, length, func) {
  if (!isFunction(func))
    func = false;

  var arr = [];
  var hival = calcRdx(min, max);
  var step = (max * hival - min * hival) / ((length - 1) * hival);
  var current = min;
  var cnt;

  // Current is assigned using a technique to compensate for IEEE error.
  // TODO: Needs better implementation.
  for (cnt = 0;
       current <= max && cnt < length;
       cnt++, current = (min * hival + step * hival * cnt) / hival) {
    arr.push((func ? func(current, cnt) : current));
  }

  return arr;
};


// arange(5) -> [0,1,2,3,4]
// arange(1,5) -> [1,2,3,4]
// arange(5,1,-1) -> [5,4,3,2]
jStat.arange = function arange(start, end, step) {
  var rl = [];
  var i;
  step = step || 1;
  if (end === undefined) {
    end = start;
    start = 0;
  }
  if (start === end || step === 0) {
    return [];
  }
  if (start < end && step < 0) {
    return [];
  }
  if (start > end && step > 0) {
    return [];
  }
  if (step > 0) {
    for (i = start; i < end; i += step) {
      rl.push(i);
    }
  } else {
    for (i = start; i > end; i += step) {
      rl.push(i);
    }
  }
  return rl;
};


// A=[[1,2,3],[4,5,6],[7,8,9]]
// slice(A,{row:{end:2},col:{start:1}}) -> [[2,3],[5,6]]
// slice(A,1,{start:1}) -> [5,6]
// as numpy code A[:2,1:]
jStat.slice = (function(){
  function _slice(list, start, end, step) {
    // note it's not equal to range.map mode it's a bug
    var i;
    var rl = [];
    var length = list.length;
    if (start === undefined && end === undefined && step === undefined) {
      return jStat.copy(list);
    }

    start = start || 0;
    end = end || list.length;
    start = start >= 0 ? start : length + start;
    end = end >= 0 ? end : length + end;
    step = step || 1;
    if (start === end || step === 0) {
      return [];
    }
    if (start < end && step < 0) {
      return [];
    }
    if (start > end && step > 0) {
      return [];
    }
    if (step > 0) {
      for (i = start; i < end; i += step) {
        rl.push(list[i]);
      }
    } else {
      for (i = start; i > end;i += step) {
        rl.push(list[i]);
      }
    }
    return rl;
  }

  function slice(list, rcSlice) {
    var colSlice, rowSlice;
    rcSlice = rcSlice || {};
    if (isNumber(rcSlice.row)) {
      if (isNumber(rcSlice.col))
        return list[rcSlice.row][rcSlice.col];
      var row = jStat.rowa(list, rcSlice.row);
      colSlice = rcSlice.col || {};
      return _slice(row, colSlice.start, colSlice.end, colSlice.step);
    }

    if (isNumber(rcSlice.col)) {
      var col = jStat.cola(list, rcSlice.col);
      rowSlice = rcSlice.row || {};
      return _slice(col, rowSlice.start, rowSlice.end, rowSlice.step);
    }

    rowSlice = rcSlice.row || {};
    colSlice = rcSlice.col || {};
    var rows = _slice(list, rowSlice.start, rowSlice.end, rowSlice.step);
    return rows.map(function(row) {
      return _slice(row, colSlice.start, colSlice.end, colSlice.step);
    });
  }

  return slice;
}());


// A=[[1,2,3],[4,5,6],[7,8,9]]
// sliceAssign(A,{row:{start:1},col:{start:1}},[[0,0],[0,0]])
// A=[[1,2,3],[4,0,0],[7,0,0]]
jStat.sliceAssign = function sliceAssign(A, rcSlice, B) {
  var nl, ml;
  if (isNumber(rcSlice.row)) {
    if (isNumber(rcSlice.col))
      return A[rcSlice.row][rcSlice.col] = B;
    rcSlice.col = rcSlice.col || {};
    rcSlice.col.start = rcSlice.col.start || 0;
    rcSlice.col.end = rcSlice.col.end || A[0].length;
    rcSlice.col.step = rcSlice.col.step || 1;
    nl = jStat.arange(rcSlice.col.start,
                          Math.min(A.length, rcSlice.col.end),
                          rcSlice.col.step);
    var m = rcSlice.row;
    nl.forEach(function(n, i) {
      A[m][n] = B[i];
    });
    return A;
  }

  if (isNumber(rcSlice.col)) {
    rcSlice.row = rcSlice.row || {};
    rcSlice.row.start = rcSlice.row.start || 0;
    rcSlice.row.end = rcSlice.row.end || A.length;
    rcSlice.row.step = rcSlice.row.step || 1;
    ml = jStat.arange(rcSlice.row.start,
                          Math.min(A[0].length, rcSlice.row.end),
                          rcSlice.row.step);
    var n = rcSlice.col;
    ml.forEach(function(m, j) {
      A[m][n] = B[j];
    });
    return A;
  }

  if (B[0].length === undefined) {
    B = [B];
  }
  rcSlice.row.start = rcSlice.row.start || 0;
  rcSlice.row.end = rcSlice.row.end || A.length;
  rcSlice.row.step = rcSlice.row.step || 1;
  rcSlice.col.start = rcSlice.col.start || 0;
  rcSlice.col.end = rcSlice.col.end || A[0].length;
  rcSlice.col.step = rcSlice.col.step || 1;
  ml = jStat.arange(rcSlice.row.start,
                        Math.min(A.length, rcSlice.row.end),
                        rcSlice.row.step);
  nl = jStat.arange(rcSlice.col.start,
                        Math.min(A[0].length, rcSlice.col.end),
                        rcSlice.col.step);
  ml.forEach(function(m, i) {
    nl.forEach(function(n, j) {
      A[m][n] = B[i][j];
    });
  });
  return A;
};


// [1,2,3] ->
// [[1,0,0],[0,2,0],[0,0,3]]
jStat.diagonal = function diagonal(diagArray) {
  var mat = jStat.zeros(diagArray.length, diagArray.length);
  diagArray.forEach(function(t, i) {
    mat[i][i] = t;
  });
  return mat;
};


// return copy of A
jStat.copy = function copy(A) {
  return A.map(function(row) {
    if (isNumber(row))
      return row;
    return row.map(function(t) {
      return t;
    });
  });
};


// TODO: Go over this entire implementation. Seems a tragic waste of resources
// doing all this work. Instead, and while ugly, use new Function() to generate
// a custom function for each static method.

// Quick reference.
var jProto = jStat.prototype;

// Default length.
jProto.length = 0;

// For internal use only.
// TODO: Check if they're actually used, and if they are then rename them
// to _*
jProto.push = Array.prototype.push;
jProto.sort = Array.prototype.sort;
jProto.splice = Array.prototype.splice;
jProto.slice = Array.prototype.slice;


// Return a clean array.
jProto.toArray = function toArray() {
  return this.length > 1 ? slice.call(this) : slice.call(this)[0];
};


// Map a function to a matrix or vector.
jProto.map = function map(func, toAlter) {
  return jStat(jStat.map(this, func, toAlter));
};


// Cumulatively combine the elements of a matrix or vector using a function.
jProto.cumreduce = function cumreduce(func, toAlter) {
  return jStat(jStat.cumreduce(this, func, toAlter));
};


// Destructively alter an array.
jProto.alter = function alter(func) {
  jStat.alter(this, func);
  return this;
};


// Extend prototype with methods that have no argument.
(function(funcs) {
  for (var i = 0; i < funcs.length; i++) (function(passfunc) {
    jProto[passfunc] = function(func) {
      var self = this,
      results;
      // Check for callback.
      if (func) {
        setTimeout(function() {
          func.call(self, jProto[passfunc].call(self));
        });
        return this;
      }
      results = jStat[passfunc](this);
      return isArray(results) ? jStat(results) : results;
    };
  })(funcs[i]);
})('transpose clear symmetric rows cols dimensions diag antidiag'.split(' '));


// Extend prototype with methods that have one argument.
(function(funcs) {
  for (var i = 0; i < funcs.length; i++) (function(passfunc) {
    jProto[passfunc] = function(index, func) {
      var self = this;
      // check for callback
      if (func) {
        setTimeout(function() {
          func.call(self, jProto[passfunc].call(self, index));
        });
        return this;
      }
      return jStat(jStat[passfunc](this, index));
    };
  })(funcs[i]);
})('row col'.split(' '));


// Extend prototype with simple shortcut methods.
(function(funcs) {
  for (var i = 0; i < funcs.length; i++) (function(passfunc) {
    jProto[passfunc] = function() {
      return jStat(jStat[passfunc].apply(null, arguments));
    };
  })(funcs[i]);
})('create zeros ones rand identity'.split(' '));


// Exposing jStat.
return jStat;

}(Math));
(function(jStat, Math) {

var isFunction = jStat.utils.isFunction;

// Ascending functions for sort
function ascNum(a, b) { return a - b; }

function clip(arg, min, max) {
  return Math.max(min, Math.min(arg, max));
}


// sum of an array
jStat.sum = function sum(arr) {
  var sum = 0;
  var i = arr.length;
  while (--i >= 0)
    sum += arr[i];
  return sum;
};


// sum squared
jStat.sumsqrd = function sumsqrd(arr) {
  var sum = 0;
  var i = arr.length;
  while (--i >= 0)
    sum += arr[i] * arr[i];
  return sum;
};


// sum of squared errors of prediction (SSE)
jStat.sumsqerr = function sumsqerr(arr) {
  var mean = jStat.mean(arr);
  var sum = 0;
  var i = arr.length;
  var tmp;
  while (--i >= 0) {
    tmp = arr[i] - mean;
    sum += tmp * tmp;
  }
  return sum;
};

// sum of an array in each row
jStat.sumrow = function sumrow(arr) {
  var sum = 0;
  var i = arr.length;
  while (--i >= 0)
    sum += arr[i];
  return sum;
};

// product of an array
jStat.product = function product(arr) {
  var prod = 1;
  var i = arr.length;
  while (--i >= 0)
    prod *= arr[i];
  return prod;
};


// minimum value of an array
jStat.min = function min(arr) {
  var low = arr[0];
  var i = 0;
  while (++i < arr.length)
    if (arr[i] < low)
      low = arr[i];
  return low;
};


// maximum value of an array
jStat.max = function max(arr) {
  var high = arr[0];
  var i = 0;
  while (++i < arr.length)
    if (arr[i] > high)
      high = arr[i];
  return high;
};


// unique values of an array
jStat.unique = function unique(arr) {
  var hash = {}, _arr = [];
  for(var i = 0; i < arr.length; i++) {
    if (!hash[arr[i]]) {
      hash[arr[i]] = true;
      _arr.push(arr[i]);
    }
  }
  return _arr;
};


// mean value of an array
jStat.mean = function mean(arr) {
  return jStat.sum(arr) / arr.length;
};


// mean squared error (MSE)
jStat.meansqerr = function meansqerr(arr) {
  return jStat.sumsqerr(arr) / arr.length;
};


// geometric mean of an array
jStat.geomean = function geomean(arr) {
  var logs = arr.map(Math.log)
  var meanOfLogs = jStat.mean(logs)
  return Math.exp(meanOfLogs)
};


// median of an array
jStat.median = function median(arr) {
  var arrlen = arr.length;
  var _arr = arr.slice().sort(ascNum);
  // check if array is even or odd, then return the appropriate
  return !(arrlen & 1)
    ? (_arr[(arrlen / 2) - 1 ] + _arr[(arrlen / 2)]) / 2
    : _arr[(arrlen / 2) | 0 ];
};


// cumulative sum of an array
jStat.cumsum = function cumsum(arr) {
  return jStat.cumreduce(arr, function (a, b) { return a + b; });
};


// cumulative product of an array
jStat.cumprod = function cumprod(arr) {
  return jStat.cumreduce(arr, function (a, b) { return a * b; });
};


// successive differences of a sequence
jStat.diff = function diff(arr) {
  var diffs = [];
  var arrLen = arr.length;
  var i;
  for (i = 1; i < arrLen; i++)
    diffs.push(arr[i] - arr[i - 1]);
  return diffs;
};


// ranks of an array
jStat.rank = function (arr) {
  var i;
  var distinctNumbers = [];
  var numberCounts = {};
  for (i = 0; i < arr.length; i++) {
    var number = arr[i];
    if (numberCounts[number]) {
      numberCounts[number]++;
    } else {
      numberCounts[number] = 1;
      distinctNumbers.push(number);
    }
  }

  var sortedDistinctNumbers = distinctNumbers.sort(ascNum);
  var numberRanks = {};
  var currentRank = 1;
  for (i = 0; i < sortedDistinctNumbers.length; i++) {
    var number = sortedDistinctNumbers[i];
    var count = numberCounts[number];
    var first = currentRank;
    var last = currentRank + count - 1;
    var rank = (first + last) / 2;
    numberRanks[number] = rank;
    currentRank += count;
  }

  return arr.map(function (number) {
    return numberRanks[number];
  });
};


// mode of an array
// if there are multiple modes of an array, return all of them
// is this the appropriate way of handling it?
jStat.mode = function mode(arr) {
  var arrLen = arr.length;
  var _arr = arr.slice().sort(ascNum);
  var count = 1;
  var maxCount = 0;
  var numMaxCount = 0;
  var mode_arr = [];
  var i;

  for (i = 0; i < arrLen; i++) {
    if (_arr[i] === _arr[i + 1]) {
      count++;
    } else {
      if (count > maxCount) {
        mode_arr = [_arr[i]];
        maxCount = count;
        numMaxCount = 0;
      }
      // are there multiple max counts
      else if (count === maxCount) {
        mode_arr.push(_arr[i]);
        numMaxCount++;
      }
      // resetting count for new value in array
      count = 1;
    }
  }

  return numMaxCount === 0 ? mode_arr[0] : mode_arr;
};


// range of an array
jStat.range = function range(arr) {
  return jStat.max(arr) - jStat.min(arr);
};

// variance of an array
// flag = true indicates sample instead of population
jStat.variance = function variance(arr, flag) {
  return jStat.sumsqerr(arr) / (arr.length - (flag ? 1 : 0));
};

// pooled variance of an array of arrays
jStat.pooledvariance = function pooledvariance(arr) {
  var sumsqerr = arr.reduce(function (a, samples) {return a + jStat.sumsqerr(samples);}, 0);
  var count = arr.reduce(function (a, samples) {return a + samples.length;}, 0);
  return sumsqerr / (count - arr.length);
};

// deviation of an array
jStat.deviation = function (arr) {
  var mean = jStat.mean(arr);
  var arrlen = arr.length;
  var dev = new Array(arrlen);
  for (var i = 0; i < arrlen; i++) {
    dev[i] = arr[i] - mean;
  }
  return dev;
};

// standard deviation of an array
// flag = true indicates sample instead of population
jStat.stdev = function stdev(arr, flag) {
  return Math.sqrt(jStat.variance(arr, flag));
};

// pooled standard deviation of an array of arrays
jStat.pooledstdev = function pooledstdev(arr) {
  return Math.sqrt(jStat.pooledvariance(arr));
};

// mean deviation (mean absolute deviation) of an array
jStat.meandev = function meandev(arr) {
  var mean = jStat.mean(arr);
  var a = [];
  for (var i = arr.length - 1; i >= 0; i--) {
    a.push(Math.abs(arr[i] - mean));
  }
  return jStat.mean(a);
};


// median deviation (median absolute deviation) of an array
jStat.meddev = function meddev(arr) {
  var median = jStat.median(arr);
  var a = [];
  for (var i = arr.length - 1; i >= 0; i--) {
    a.push(Math.abs(arr[i] - median));
  }
  return jStat.median(a);
};


// coefficient of variation
jStat.coeffvar = function coeffvar(arr) {
  return jStat.stdev(arr) / jStat.mean(arr);
};


// quartiles of an array
jStat.quartiles = function quartiles(arr) {
  var arrlen = arr.length;
  var _arr = arr.slice().sort(ascNum);
  return [
    _arr[ Math.round((arrlen) / 4) - 1 ],
    _arr[ Math.round((arrlen) / 2) - 1 ],
    _arr[ Math.round((arrlen) * 3 / 4) - 1 ]
  ];
};


// Arbitary quantiles of an array. Direct port of the scipy.stats
// implementation by Pierre GF Gerard-Marchant.
jStat.quantiles = function quantiles(arr, quantilesArray, alphap, betap) {
  var sortedArray = arr.slice().sort(ascNum);
  var quantileVals = [quantilesArray.length];
  var n = arr.length;
  var i, p, m, aleph, k, gamma;

  if (typeof alphap === 'undefined')
    alphap = 3 / 8;
  if (typeof betap === 'undefined')
    betap = 3 / 8;

  for (i = 0; i < quantilesArray.length; i++) {
    p = quantilesArray[i];
    m = alphap + p * (1 - alphap - betap);
    aleph = n * p + m;
    k = Math.floor(clip(aleph, 1, n - 1));
    gamma = clip(aleph - k, 0, 1);
    quantileVals[i] = (1 - gamma) * sortedArray[k - 1] + gamma * sortedArray[k];
  }

  return quantileVals;
};

// Return the k-th percentile of values in a range, where k is in the range 0..1, inclusive.
// Passing true for the exclusive parameter excludes both endpoints of the range.
jStat.percentile = function percentile(arr, k, exclusive) {
  var _arr = arr.slice().sort(ascNum);
  var realIndex = k * (_arr.length + (exclusive ? 1 : -1)) + (exclusive ? 0 : 1);
  var index = parseInt(realIndex);
  var frac = realIndex - index;
  if (index + 1 < _arr.length) {
    return _arr[index - 1] + frac * (_arr[index] - _arr[index - 1]);
  } else {
    return _arr[index - 1];
  }
}

// The percentile rank of score in a given array. Returns the percentage
// of all values in the input array that are less than (kind='strict') or
// less or equal than (kind='weak') score. Default is weak.
jStat.percentileOfScore = function percentileOfScore(arr, score, kind) {
  var counter = 0;
  var len = arr.length;
  var strict = false;
  var value, i;

  if (kind === 'strict')
    strict = true;

  for (i = 0; i < len; i++) {
    value = arr[i];
    if ((strict && value < score) ||
        (!strict && value <= score)) {
      counter++;
    }
  }

  return counter / len;
};


// Histogram (bin count) data
jStat.histogram = function histogram(arr, binCnt) {
  binCnt = binCnt || 4;
  var first = jStat.min(arr);
  var binWidth = (jStat.max(arr) - first) / binCnt;
  var len = arr.length;
  var bins = [];
  var i;

  for (i = 0; i < binCnt; i++)
    bins[i] = 0;
  for (i = 0; i < len; i++)
    bins[Math.min(Math.floor(((arr[i] - first) / binWidth)), binCnt - 1)] += 1;

  return bins;
};


// covariance of two arrays
jStat.covariance = function covariance(arr1, arr2) {
  var u = jStat.mean(arr1);
  var v = jStat.mean(arr2);
  var arr1Len = arr1.length;
  var sq_dev = new Array(arr1Len);
  var i;

  for (i = 0; i < arr1Len; i++)
    sq_dev[i] = (arr1[i] - u) * (arr2[i] - v);

  return jStat.sum(sq_dev) / (arr1Len - 1);
};


// (pearson's) population correlation coefficient, rho
jStat.corrcoeff = function corrcoeff(arr1, arr2) {
  return jStat.covariance(arr1, arr2) /
      jStat.stdev(arr1, 1) /
      jStat.stdev(arr2, 1);
};

  // (spearman's) rank correlation coefficient, sp
jStat.spearmancoeff =  function (arr1, arr2) {
  arr1 = jStat.rank(arr1);
  arr2 = jStat.rank(arr2);
  //return pearson's correlation of the ranks:
  return jStat.corrcoeff(arr1, arr2);
}


// statistical standardized moments (general form of skew/kurt)
jStat.stanMoment = function stanMoment(arr, n) {
  var mu = jStat.mean(arr);
  var sigma = jStat.stdev(arr);
  var len = arr.length;
  var skewSum = 0;

  for (var i = 0; i < len; i++)
    skewSum += Math.pow((arr[i] - mu) / sigma, n);

  return skewSum / arr.length;
};

// (pearson's) moment coefficient of skewness
jStat.skewness = function skewness(arr) {
  return jStat.stanMoment(arr, 3);
};

// (pearson's) (excess) kurtosis
jStat.kurtosis = function kurtosis(arr) {
  return jStat.stanMoment(arr, 4) - 3;
};


var jProto = jStat.prototype;


// Extend jProto with method for calculating cumulative sums and products.
// This differs from the similar extension below as cumsum and cumprod should
// not be run again in the case fullbool === true.
// If a matrix is passed, automatically assume operation should be done on the
// columns.
(function(funcs) {
  for (var i = 0; i < funcs.length; i++) (function(passfunc) {
    // If a matrix is passed, automatically assume operation should be done on
    // the columns.
    jProto[passfunc] = function(fullbool, func) {
      var arr = [];
      var i = 0;
      var tmpthis = this;
      // Assignment reassignation depending on how parameters were passed in.
      if (isFunction(fullbool)) {
        func = fullbool;
        fullbool = false;
      }
      // Check if a callback was passed with the function.
      if (func) {
        setTimeout(function() {
          func.call(tmpthis, jProto[passfunc].call(tmpthis, fullbool));
        });
        return this;
      }
      // Check if matrix and run calculations.
      if (this.length > 1) {
        tmpthis = fullbool === true ? this : this.transpose();
        for (; i < tmpthis.length; i++)
          arr[i] = jStat[passfunc](tmpthis[i]);
        return arr;
      }
      // Pass fullbool if only vector, not a matrix. for variance and stdev.
      return jStat[passfunc](this[0], fullbool);
    };
  })(funcs[i]);
})(('cumsum cumprod').split(' '));


// Extend jProto with methods which don't require arguments and work on columns.
(function(funcs) {
  for (var i = 0; i < funcs.length; i++) (function(passfunc) {
    // If a matrix is passed, automatically assume operation should be done on
    // the columns.
    jProto[passfunc] = function(fullbool, func) {
      var arr = [];
      var i = 0;
      var tmpthis = this;
      // Assignment reassignation depending on how parameters were passed in.
      if (isFunction(fullbool)) {
        func = fullbool;
        fullbool = false;
      }
      // Check if a callback was passed with the function.
      if (func) {
        setTimeout(function() {
          func.call(tmpthis, jProto[passfunc].call(tmpthis, fullbool));
        });
        return this;
      }
      // Check if matrix and run calculations.
      if (this.length > 1) {
        if (passfunc !== 'sumrow')
          tmpthis = fullbool === true ? this : this.transpose();
        for (; i < tmpthis.length; i++)
          arr[i] = jStat[passfunc](tmpthis[i]);
        return fullbool === true
            ? jStat[passfunc](jStat.utils.toVector(arr))
            : arr;
      }
      // Pass fullbool if only vector, not a matrix. for variance and stdev.
      return jStat[passfunc](this[0], fullbool);
    };
  })(funcs[i]);
})(('sum sumsqrd sumsqerr sumrow product min max unique mean meansqerr ' +
    'geomean median diff rank mode range variance deviation stdev meandev ' +
    'meddev coeffvar quartiles histogram skewness kurtosis').split(' '));


// Extend jProto with functions that take arguments. Operations on matrices are
// done on columns.
(function(funcs) {
  for (var i = 0; i < funcs.length; i++) (function(passfunc) {
    jProto[passfunc] = function() {
      var arr = [];
      var i = 0;
      var tmpthis = this;
      var args = Array.prototype.slice.call(arguments);
      var callbackFunction;

      // If the last argument is a function, we assume it's a callback; we
      // strip the callback out and call the function again.
      if (isFunction(args[args.length - 1])) {
        callbackFunction = args[args.length - 1];
        var argsToPass = args.slice(0, args.length - 1);

        setTimeout(function() {
          callbackFunction.call(tmpthis,
                                jProto[passfunc].apply(tmpthis, argsToPass));
        });
        return this;

      // Otherwise we curry the function args and call normally.
      } else {
        callbackFunction = undefined;
        var curriedFunction = function curriedFunction(vector) {
          return jStat[passfunc].apply(tmpthis, [vector].concat(args));
        }
      }

      // If this is a matrix, run column-by-column.
      if (this.length > 1) {
        tmpthis = tmpthis.transpose();
        for (; i < tmpthis.length; i++)
          arr[i] = curriedFunction(tmpthis[i]);
        return arr;
      }

      // Otherwise run on the vector.
      return curriedFunction(this[0]);
    };
  })(funcs[i]);
})('quantiles percentileOfScore'.split(' '));

}(jStat, Math));
// Special functions //
(function(jStat, Math) {

// Log-gamma function
jStat.gammaln = function gammaln(x) {
  var j = 0;
  var cof = [
    76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
  ];
  var ser = 1.000000000190015;
  var xx, y, tmp;
  tmp = (y = xx = x) + 5.5;
  tmp -= (xx + 0.5) * Math.log(tmp);
  for (; j < 6; j++)
    ser += cof[j] / ++y;
  return Math.log(2.5066282746310005 * ser / xx) - tmp;
};

/*
 * log-gamma function to support poisson distribution sampling. The
 * algorithm comes from SPECFUN by Shanjie Zhang and Jianming Jin and their
 * book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
 */
jStat.loggam = function loggam(x) {
  var x0, x2, xp, gl, gl0;
  var k, n;

  var a = [8.333333333333333e-02, -2.777777777777778e-03,
          7.936507936507937e-04, -5.952380952380952e-04,
          8.417508417508418e-04, -1.917526917526918e-03,
          6.410256410256410e-03, -2.955065359477124e-02,
          1.796443723688307e-01, -1.39243221690590e+00];
  x0 = x;
  n = 0;
  if ((x == 1.0) || (x == 2.0)) {
      return 0.0;
  }
  if (x <= 7.0) {
      n = Math.floor(7 - x);
      x0 = x + n;
  }
  x2 = 1.0 / (x0 * x0);
  xp = 2 * Math.PI;
  gl0 = a[9];
  for (k = 8; k >= 0; k--) {
      gl0 *= x2;
      gl0 += a[k];
  }
  gl = gl0 / x0 + 0.5 * Math.log(xp) + (x0 - 0.5) * Math.log(x0) - x0;
  if (x <= 7.0) {
      for (k = 1; k <= n; k++) {
          gl -= Math.log(x0 - 1.0);
          x0 -= 1.0;
      }
  }
  return gl;
}

// gamma of x
jStat.gammafn = function gammafn(x) {
  var p = [-1.716185138865495, 24.76565080557592, -379.80425647094563,
           629.3311553128184, 866.9662027904133, -31451.272968848367,
           -36144.413418691176, 66456.14382024054
  ];
  var q = [-30.8402300119739, 315.35062697960416, -1015.1563674902192,
           -3107.771671572311, 22538.118420980151, 4755.8462775278811,
           -134659.9598649693, -115132.2596755535];
  var fact = false;
  var n = 0;
  var xden = 0;
  var xnum = 0;
  var y = x;
  var i, z, yi, res;
  if (x > 171.6243769536076) {
    return Infinity;
  }
  if (y <= 0) {
    res = y % 1 + 3.6e-16;
    if (res) {
      fact = (!(y & 1) ? 1 : -1) * Math.PI / Math.sin(Math.PI * res);
      y = 1 - y;
    } else {
      return Infinity;
    }
  }
  yi = y;
  if (y < 1) {
    z = y++;
  } else {
    z = (y -= n = (y | 0) - 1) - 1;
  }
  for (i = 0; i < 8; ++i) {
    xnum = (xnum + p[i]) * z;
    xden = xden * z + q[i];
  }
  res = xnum / xden + 1;
  if (yi < y) {
    res /= yi;
  } else if (yi > y) {
    for (i = 0; i < n; ++i) {
      res *= y;
      y++;
    }
  }
  if (fact) {
    res = fact / res;
  }
  return res;
};


// lower incomplete gamma function, which is usually typeset with a
// lower-case greek gamma as the function symbol
jStat.gammap = function gammap(a, x) {
  return jStat.lowRegGamma(a, x) * jStat.gammafn(a);
};


// The lower regularized incomplete gamma function, usually written P(a,x)
jStat.lowRegGamma = function lowRegGamma(a, x) {
  var aln = jStat.gammaln(a);
  var ap = a;
  var sum = 1 / a;
  var del = sum;
  var b = x + 1 - a;
  var c = 1 / 1.0e-30;
  var d = 1 / b;
  var h = d;
  var i = 1;
  // calculate maximum number of itterations required for a
  var ITMAX = -~(Math.log((a >= 1) ? a : 1 / a) * 8.5 + a * 0.4 + 17);
  var an;

  if (x < 0 || a <= 0) {
    return NaN;
  } else if (x < a + 1) {
    for (; i <= ITMAX; i++) {
      sum += del *= x / ++ap;
    }
    return (sum * Math.exp(-x + a * Math.log(x) - (aln)));
  }

  for (; i <= ITMAX; i++) {
    an = -i * (i - a);
    b += 2;
    d = an * d + b;
    c = b + an / c;
    d = 1 / d;
    h *= d * c;
  }

  return (1 - h * Math.exp(-x + a * Math.log(x) - (aln)));
};

// natural log factorial of n
jStat.factorialln = function factorialln(n) {
  return n < 0 ? NaN : jStat.gammaln(n + 1);
};

// factorial of n
jStat.factorial = function factorial(n) {
  return n < 0 ? NaN : jStat.gammafn(n + 1);
};

// combinations of n, m
jStat.combination = function combination(n, m) {
  // make sure n or m don't exceed the upper limit of usable values
  return (n > 170 || m > 170)
      ? Math.exp(jStat.combinationln(n, m))
      : (jStat.factorial(n) / jStat.factorial(m)) / jStat.factorial(n - m);
};


jStat.combinationln = function combinationln(n, m){
  return jStat.factorialln(n) - jStat.factorialln(m) - jStat.factorialln(n - m);
};


// permutations of n, m
jStat.permutation = function permutation(n, m) {
  return jStat.factorial(n) / jStat.factorial(n - m);
};


// beta function
jStat.betafn = function betafn(x, y) {
  // ensure arguments are positive
  if (x <= 0 || y <= 0)
    return undefined;
  // make sure x + y doesn't exceed the upper limit of usable values
  return (x + y > 170)
      ? Math.exp(jStat.betaln(x, y))
      : jStat.gammafn(x) * jStat.gammafn(y) / jStat.gammafn(x + y);
};


// natural logarithm of beta function
jStat.betaln = function betaln(x, y) {
  return jStat.gammaln(x) + jStat.gammaln(y) - jStat.gammaln(x + y);
};


// Evaluates the continued fraction for incomplete beta function by modified
// Lentz's method.
jStat.betacf = function betacf(x, a, b) {
  var fpmin = 1e-30;
  var m = 1;
  var qab = a + b;
  var qap = a + 1;
  var qam = a - 1;
  var c = 1;
  var d = 1 - qab * x / qap;
  var m2, aa, del, h;

  // These q's will be used in factors that occur in the coefficients
  if (Math.abs(d) < fpmin)
    d = fpmin;
  d = 1 / d;
  h = d;

  for (; m <= 100; m++) {
    m2 = 2 * m;
    aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    // One step (the even one) of the recurrence
    d = 1 + aa * d;
    if (Math.abs(d) < fpmin)
      d = fpmin;
    c = 1 + aa / c;
    if (Math.abs(c) < fpmin)
      c = fpmin;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    // Next step of the recurrence (the odd one)
    d = 1 + aa * d;
    if (Math.abs(d) < fpmin)
      d = fpmin;
    c = 1 + aa / c;
    if (Math.abs(c) < fpmin)
      c = fpmin;
    d = 1 / d;
    del = d * c;
    h *= del;
    if (Math.abs(del - 1.0) < 3e-7)
      break;
  }

  return h;
};


// Returns the inverse of the lower regularized inomplete gamma function
jStat.gammapinv = function gammapinv(p, a) {
  var j = 0;
  var a1 = a - 1;
  var EPS = 1e-8;
  var gln = jStat.gammaln(a);
  var x, err, t, u, pp, lna1, afac;

  if (p >= 1)
    return Math.max(100, a + 100 * Math.sqrt(a));
  if (p <= 0)
    return 0;
  if (a > 1) {
    lna1 = Math.log(a1);
    afac = Math.exp(a1 * (lna1 - 1) - gln);
    pp = (p < 0.5) ? p : 1 - p;
    t = Math.sqrt(-2 * Math.log(pp));
    x = (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t;
    if (p < 0.5)
      x = -x;
    x = Math.max(1e-3,
                 a * Math.pow(1 - 1 / (9 * a) - x / (3 * Math.sqrt(a)), 3));
  } else {
    t = 1 - a * (0.253 + a * 0.12);
    if (p < t)
      x = Math.pow(p / t, 1 / a);
    else
      x = 1 - Math.log(1 - (p - t) / (1 - t));
  }

  for(; j < 12; j++) {
    if (x <= 0)
      return 0;
    err = jStat.lowRegGamma(a, x) - p;
    if (a > 1)
      t = afac * Math.exp(-(x - a1) + a1 * (Math.log(x) - lna1));
    else
      t = Math.exp(-x + a1 * Math.log(x) - gln);
    u = err / t;
    x -= (t = u / (1 - 0.5 * Math.min(1, u * ((a - 1) / x - 1))));
    if (x <= 0)
      x = 0.5 * (x + t);
    if (Math.abs(t) < EPS * x)
      break;
  }

  return x;
};


// Returns the error function erf(x)
jStat.erf = function erf(x) {
  var cof = [-1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2,
             -9.561514786808631e-3, -9.46595344482036e-4, 3.66839497852761e-4,
             4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6,
             1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8,
             6.529054439e-9, 5.059343495e-9, -9.91364156e-10,
             -2.27365122e-10, 9.6467911e-11, 2.394038e-12,
             -6.886027e-12, 8.94487e-13, 3.13092e-13,
             -1.12708e-13, 3.81e-16, 7.106e-15,
             -1.523e-15, -9.4e-17, 1.21e-16,
             -2.8e-17];
  var j = cof.length - 1;
  var isneg = false;
  var d = 0;
  var dd = 0;
  var t, ty, tmp, res;

  if (x < 0) {
    x = -x;
    isneg = true;
  }

  t = 2 / (2 + x);
  ty = 4 * t - 2;

  for(; j > 0; j--) {
    tmp = d;
    d = ty * d - dd + cof[j];
    dd = tmp;
  }

  res = t * Math.exp(-x * x + 0.5 * (cof[0] + ty * d) - dd);
  return isneg ? res - 1 : 1 - res;
};


// Returns the complmentary error function erfc(x)
jStat.erfc = function erfc(x) {
  return 1 - jStat.erf(x);
};


// Returns the inverse of the complementary error function
jStat.erfcinv = function erfcinv(p) {
  var j = 0;
  var x, err, t, pp;
  if (p >= 2)
    return -100;
  if (p <= 0)
    return 100;
  pp = (p < 1) ? p : 2 - p;
  t = Math.sqrt(-2 * Math.log(pp / 2));
  x = -0.70711 * ((2.30753 + t * 0.27061) /
                  (1 + t * (0.99229 + t * 0.04481)) - t);
  for (; j < 2; j++) {
    err = jStat.erfc(x) - pp;
    x += err / (1.12837916709551257 * Math.exp(-x * x) - x * err);
  }
  return (p < 1) ? x : -x;
};


// Returns the inverse of the incomplete beta function
jStat.ibetainv = function ibetainv(p, a, b) {
  var EPS = 1e-8;
  var a1 = a - 1;
  var b1 = b - 1;
  var j = 0;
  var lna, lnb, pp, t, u, err, x, al, h, w, afac;
  if (p <= 0)
    return 0;
  if (p >= 1)
    return 1;
  if (a >= 1 && b >= 1) {
    pp = (p < 0.5) ? p : 1 - p;
    t = Math.sqrt(-2 * Math.log(pp));
    x = (2.30753 + t * 0.27061) / (1 + t* (0.99229 + t * 0.04481)) - t;
    if (p < 0.5)
      x = -x;
    al = (x * x - 3) / 6;
    h = 2 / (1 / (2 * a - 1)  + 1 / (2 * b - 1));
    w = (x * Math.sqrt(al + h) / h) - (1 / (2 * b - 1) - 1 / (2 * a - 1)) *
        (al + 5 / 6 - 2 / (3 * h));
    x = a / (a + b * Math.exp(2 * w));
  } else {
    lna = Math.log(a / (a + b));
    lnb = Math.log(b / (a + b));
    t = Math.exp(a * lna) / a;
    u = Math.exp(b * lnb) / b;
    w = t + u;
    if (p < t / w)
      x = Math.pow(a * w * p, 1 / a);
    else
      x = 1 - Math.pow(b * w * (1 - p), 1 / b);
  }
  afac = -jStat.gammaln(a) - jStat.gammaln(b) + jStat.gammaln(a + b);
  for(; j < 10; j++) {
    if (x === 0 || x === 1)
      return x;
    err = jStat.ibeta(x, a, b) - p;
    t = Math.exp(a1 * Math.log(x) + b1 * Math.log(1 - x) + afac);
    u = err / t;
    x -= (t = u / (1 - 0.5 * Math.min(1, u * (a1 / x - b1 / (1 - x)))));
    if (x <= 0)
      x = 0.5 * (x + t);
    if (x >= 1)
      x = 0.5 * (x + t + 1);
    if (Math.abs(t) < EPS * x && j > 0)
      break;
  }
  return x;
};


// Returns the incomplete beta function I_x(a,b)
jStat.ibeta = function ibeta(x, a, b) {
  // Factors in front of the continued fraction.
  var bt = (x === 0 || x === 1) ?  0 :
    Math.exp(jStat.gammaln(a + b) - jStat.gammaln(a) -
             jStat.gammaln(b) + a * Math.log(x) + b *
             Math.log(1 - x));
  if (x < 0 || x > 1)
    return false;
  if (x < (a + 1) / (a + b + 2))
    // Use continued fraction directly.
    return bt * jStat.betacf(x, a, b) / a;
  // else use continued fraction after making the symmetry transformation.
  return 1 - bt * jStat.betacf(1 - x, b, a) / b;
};


// Returns a normal deviate (mu=0, sigma=1).
// If n and m are specified it returns a object of normal deviates.
jStat.randn = function randn(n, m) {
  var u, v, x, y, q;
  if (!m)
    m = n;
  if (n)
    return jStat.create(n, m, function() { return jStat.randn(); });
  do {
    u = jStat._random_fn();
    v = 1.7156 * (jStat._random_fn() - 0.5);
    x = u - 0.449871;
    y = Math.abs(v) + 0.386595;
    q = x * x + y * (0.19600 * y - 0.25472 * x);
  } while (q > 0.27597 && (q > 0.27846 || v * v > -4 * Math.log(u) * u * u));
  return v / u;
};


// Returns a gamma deviate by the method of Marsaglia and Tsang.
jStat.randg = function randg(shape, n, m) {
  var oalph = shape;
  var a1, a2, u, v, x, mat;
  if (!m)
    m = n;
  if (!shape)
    shape = 1;
  if (n) {
    mat = jStat.zeros(n,m);
    mat.alter(function() { return jStat.randg(shape); });
    return mat;
  }
  if (shape < 1)
    shape += 1;
  a1 = shape - 1 / 3;
  a2 = 1 / Math.sqrt(9 * a1);
  do {
    do {
      x = jStat.randn();
      v = 1 + a2 * x;
    } while(v <= 0);
    v = v * v * v;
    u = jStat._random_fn();
  } while(u > 1 - 0.331 * Math.pow(x, 4) &&
          Math.log(u) > 0.5 * x*x + a1 * (1 - v + Math.log(v)));
  // alpha > 1
  if (shape == oalph)
    return a1 * v;
  // alpha < 1
  do {
    u = jStat._random_fn();
  } while(u === 0);
  return Math.pow(u, 1 / oalph) * a1 * v;
};


// making use of static methods on the instance
(function(funcs) {
  for (var i = 0; i < funcs.length; i++) (function(passfunc) {
    jStat.fn[passfunc] = function() {
      return jStat(
          jStat.map(this, function(value) { return jStat[passfunc](value); }));
    }
  })(funcs[i]);
})('gammaln gammafn factorial factorialln'.split(' '));


(function(funcs) {
  for (var i = 0; i < funcs.length; i++) (function(passfunc) {
    jStat.fn[passfunc] = function() {
      return jStat(jStat[passfunc].apply(null, arguments));
    };
  })(funcs[i]);
})('randn'.split(' '));

}(jStat, Math));
(function(jStat, Math) {

// generate all distribution instance methods
(function(list) {
  for (var i = 0; i < list.length; i++) (function(func) {
    // distribution instance method
    jStat[func] = function f(a, b, c) {
      if (!(this instanceof f))
        return new f(a, b, c);
      this._a = a;
      this._b = b;
      this._c = c;
      return this;
    };
    // distribution method to be used on a jStat instance
    jStat.fn[func] = function(a, b, c) {
      var newthis = jStat[func](a, b, c);
      newthis.data = this;
      return newthis;
    };
    // sample instance method
    jStat[func].prototype.sample = function(arr) {
      var a = this._a;
      var b = this._b;
      var c = this._c;
      if (arr)
        return jStat.alter(arr, function() {
          return jStat[func].sample(a, b, c);
        });
      else
        return jStat[func].sample(a, b, c);
    };
    // generate the pdf, cdf and inv instance methods
    (function(vals) {
      for (var i = 0; i < vals.length; i++) (function(fnfunc) {
        jStat[func].prototype[fnfunc] = function(x) {
          var a = this._a;
          var b = this._b;
          var c = this._c;
          if (!x && x !== 0)
            x = this.data;
          if (typeof x !== 'number') {
            return jStat.fn.map.call(x, function(x) {
              return jStat[func][fnfunc](x, a, b, c);
            });
          }
          return jStat[func][fnfunc](x, a, b, c);
        };
      })(vals[i]);
    })('pdf cdf inv'.split(' '));
    // generate the mean, median, mode and variance instance methods
    (function(vals) {
      for (var i = 0; i < vals.length; i++) (function(fnfunc) {
        jStat[func].prototype[fnfunc] = function() {
          return jStat[func][fnfunc](this._a, this._b, this._c);
        };
      })(vals[i]);
    })('mean median mode variance'.split(' '));
  })(list[i]);
})((
  'beta centralF cauchy chisquare exponential gamma invgamma kumaraswamy ' +
  'laplace lognormal noncentralt normal pareto studentt weibull uniform ' +
  'binomial negbin hypgeom poisson triangular tukey arcsine'
).split(' '));



// extend beta function with static methods
jStat.extend(jStat.beta, {
  pdf: function pdf(x, alpha, beta) {
    // PDF is zero outside the support
    if (x > 1 || x < 0)
      return 0;
    // PDF is one for the uniform case
    if (alpha == 1 && beta == 1)
      return 1;

    if (alpha < 512 && beta < 512) {
      return (Math.pow(x, alpha - 1) * Math.pow(1 - x, beta - 1)) /
          jStat.betafn(alpha, beta);
    } else {
      return Math.exp((alpha - 1) * Math.log(x) +
                      (beta - 1) * Math.log(1 - x) -
                      jStat.betaln(alpha, beta));
    }
  },

  cdf: function cdf(x, alpha, beta) {
    return (x > 1 || x < 0) ? (x > 1) * 1 : jStat.ibeta(x, alpha, beta);
  },

  inv: function inv(x, alpha, beta) {
    return jStat.ibetainv(x, alpha, beta);
  },

  mean: function mean(alpha, beta) {
    return alpha / (alpha + beta);
  },

  median: function median(alpha, beta) {
    return jStat.ibetainv(0.5, alpha, beta);
  },

  mode: function mode(alpha, beta) {
    return (alpha - 1 ) / ( alpha + beta - 2);
  },

  // return a random sample
  sample: function sample(alpha, beta) {
    var u = jStat.randg(alpha);
    return u / (u + jStat.randg(beta));
  },

  variance: function variance(alpha, beta) {
    return (alpha * beta) / (Math.pow(alpha + beta, 2) * (alpha + beta + 1));
  }
});

// extend F function with static methods
jStat.extend(jStat.centralF, {
  // This implementation of the pdf function avoids float overflow
  // See the way that R calculates this value:
  // https://svn.r-project.org/R/trunk/src/nmath/df.c
  pdf: function pdf(x, df1, df2) {
    var p, q, f;

    if (x < 0)
      return 0;

    if (df1 <= 2) {
      if (x === 0 && df1 < 2) {
        return Infinity;
      }
      if (x === 0 && df1 === 2) {
        return 1;
      }
      return (1 / jStat.betafn(df1 / 2, df2 / 2)) *
              Math.pow(df1 / df2, df1 / 2) *
              Math.pow(x, (df1/2) - 1) *
              Math.pow((1 + (df1 / df2) * x), -(df1 + df2) / 2);
    }

    p = (df1 * x) / (df2 + x * df1);
    q = df2 / (df2 + x * df1);
    f = df1 * q / 2.0;
    return f * jStat.binomial.pdf((df1 - 2) / 2, (df1 + df2 - 2) / 2, p);
  },

  cdf: function cdf(x, df1, df2) {
    if (x < 0)
      return 0;
    return jStat.ibeta((df1 * x) / (df1 * x + df2), df1 / 2, df2 / 2);
  },

  inv: function inv(x, df1, df2) {
    return df2 / (df1 * (1 / jStat.ibetainv(x, df1 / 2, df2 / 2) - 1));
  },

  mean: function mean(df1, df2) {
    return (df2 > 2) ? df2 / (df2 - 2) : undefined;
  },

  mode: function mode(df1, df2) {
    return (df1 > 2) ? (df2 * (df1 - 2)) / (df1 * (df2 + 2)) : undefined;
  },

  // return a random sample
  sample: function sample(df1, df2) {
    var x1 = jStat.randg(df1 / 2) * 2;
    var x2 = jStat.randg(df2 / 2) * 2;
    return (x1 / df1) / (x2 / df2);
  },

  variance: function variance(df1, df2) {
    if (df2 <= 4)
      return undefined;
    return 2 * df2 * df2 * (df1 + df2 - 2) /
        (df1 * (df2 - 2) * (df2 - 2) * (df2 - 4));
  }
});


// extend cauchy function with static methods
jStat.extend(jStat.cauchy, {
  pdf: function pdf(x, local, scale) {
    if (scale < 0) { return 0; }

    return (scale / (Math.pow(x - local, 2) + Math.pow(scale, 2))) / Math.PI;
  },

  cdf: function cdf(x, local, scale) {
    return Math.atan((x - local) / scale) / Math.PI + 0.5;
  },

  inv: function(p, local, scale) {
    return local + scale * Math.tan(Math.PI * (p - 0.5));
  },

  median: function median(local/*, scale*/) {
    return local;
  },

  mode: function mode(local/*, scale*/) {
    return local;
  },

  sample: function sample(local, scale) {
    return jStat.randn() *
        Math.sqrt(1 / (2 * jStat.randg(0.5))) * scale + local;
  }
});



// extend chisquare function with static methods
jStat.extend(jStat.chisquare, {
  pdf: function pdf(x, dof) {
    if (x < 0)
      return 0;
    return (x === 0 && dof === 2) ? 0.5 :
        Math.exp((dof / 2 - 1) * Math.log(x) - x / 2 - (dof / 2) *
                 Math.log(2) - jStat.gammaln(dof / 2));
  },

  cdf: function cdf(x, dof) {
    if (x < 0)
      return 0;
    return jStat.lowRegGamma(dof / 2, x / 2);
  },

  inv: function(p, dof) {
    return 2 * jStat.gammapinv(p, 0.5 * dof);
  },

  mean : function(dof) {
    return dof;
  },

  // TODO: this is an approximation (is there a better way?)
  median: function median(dof) {
    return dof * Math.pow(1 - (2 / (9 * dof)), 3);
  },

  mode: function mode(dof) {
    return (dof - 2 > 0) ? dof - 2 : 0;
  },

  sample: function sample(dof) {
    return jStat.randg(dof / 2) * 2;
  },

  variance: function variance(dof) {
    return 2 * dof;
  }
});



// extend exponential function with static methods
jStat.extend(jStat.exponential, {
  pdf: function pdf(x, rate) {
    return x < 0 ? 0 : rate * Math.exp(-rate * x);
  },

  cdf: function cdf(x, rate) {
    return x < 0 ? 0 : 1 - Math.exp(-rate * x);
  },

  inv: function(p, rate) {
    return -Math.log(1 - p) / rate;
  },

  mean : function(rate) {
    return 1 / rate;
  },

  median: function (rate) {
    return (1 / rate) * Math.log(2);
  },

  mode: function mode(/*rate*/) {
    return 0;
  },

  sample: function sample(rate) {
    return -1 / rate * Math.log(jStat._random_fn());
  },

  variance : function(rate) {
    return Math.pow(rate, -2);
  }
});



// extend gamma function with static methods
jStat.extend(jStat.gamma, {
  pdf: function pdf(x, shape, scale) {
    if (x < 0)
      return 0;
    return (x === 0 && shape === 1) ? 1 / scale :
            Math.exp((shape - 1) * Math.log(x) - x / scale -
                    jStat.gammaln(shape) - shape * Math.log(scale));
  },

  cdf: function cdf(x, shape, scale) {
    if (x < 0)
      return 0;
    return jStat.lowRegGamma(shape, x / scale);
  },

  inv: function(p, shape, scale) {
    return jStat.gammapinv(p, shape) * scale;
  },

  mean : function(shape, scale) {
    return shape * scale;
  },

  mode: function mode(shape, scale) {
    if(shape > 1) return (shape - 1) * scale;
    return undefined;
  },

  sample: function sample(shape, scale) {
    return jStat.randg(shape) * scale;
  },

  variance: function variance(shape, scale) {
    return shape * scale * scale;
  }
});

// extend inverse gamma function with static methods
jStat.extend(jStat.invgamma, {
  pdf: function pdf(x, shape, scale) {
    if (x <= 0)
      return 0;
    return Math.exp(-(shape + 1) * Math.log(x) - scale / x -
                    jStat.gammaln(shape) + shape * Math.log(scale));
  },

  cdf: function cdf(x, shape, scale) {
    if (x <= 0)
      return 0;
    return 1 - jStat.lowRegGamma(shape, scale / x);
  },

  inv: function(p, shape, scale) {
    return scale / jStat.gammapinv(1 - p, shape);
  },

  mean : function(shape, scale) {
    return (shape > 1) ? scale / (shape - 1) : undefined;
  },

  mode: function mode(shape, scale) {
    return scale / (shape + 1);
  },

  sample: function sample(shape, scale) {
    return scale / jStat.randg(shape);
  },

  variance: function variance(shape, scale) {
    if (shape <= 2)
      return undefined;
    return scale * scale / ((shape - 1) * (shape - 1) * (shape - 2));
  }
});


// extend kumaraswamy function with static methods
jStat.extend(jStat.kumaraswamy, {
  pdf: function pdf(x, alpha, beta) {
    if (x === 0 && alpha === 1)
      return beta;
    else if (x === 1 && beta === 1)
      return alpha;
    return Math.exp(Math.log(alpha) + Math.log(beta) + (alpha - 1) *
                    Math.log(x) + (beta - 1) *
                    Math.log(1 - Math.pow(x, alpha)));
  },

  cdf: function cdf(x, alpha, beta) {
    if (x < 0)
      return 0;
    else if (x > 1)
      return 1;
    return (1 - Math.pow(1 - Math.pow(x, alpha), beta));
  },

  inv: function inv(p, alpha, beta) {
    return Math.pow(1 - Math.pow(1 - p, 1 / beta), 1 / alpha);
  },

  mean : function(alpha, beta) {
    return (beta * jStat.gammafn(1 + 1 / alpha) *
            jStat.gammafn(beta)) / (jStat.gammafn(1 + 1 / alpha + beta));
  },

  median: function median(alpha, beta) {
    return Math.pow(1 - Math.pow(2, -1 / beta), 1 / alpha);
  },

  mode: function mode(alpha, beta) {
    if (!(alpha >= 1 && beta >= 1 && (alpha !== 1 && beta !== 1)))
      return undefined;
    return Math.pow((alpha - 1) / (alpha * beta - 1), 1 / alpha);
  },

  variance: function variance(/*alpha, beta*/) {
    throw new Error('variance not yet implemented');
    // TODO: complete this
  }
});



// extend lognormal function with static methods
jStat.extend(jStat.lognormal, {
  pdf: function pdf(x, mu, sigma) {
    if (x <= 0)
      return 0;
    return Math.exp(-Math.log(x) - 0.5 * Math.log(2 * Math.PI) -
                    Math.log(sigma) - Math.pow(Math.log(x) - mu, 2) /
                    (2 * sigma * sigma));
  },

  cdf: function cdf(x, mu, sigma) {
    if (x < 0)
      return 0;
    return 0.5 +
        (0.5 * jStat.erf((Math.log(x) - mu) / Math.sqrt(2 * sigma * sigma)));
  },

  inv: function(p, mu, sigma) {
    return Math.exp(-1.41421356237309505 * sigma * jStat.erfcinv(2 * p) + mu);
  },

  mean: function mean(mu, sigma) {
    return Math.exp(mu + sigma * sigma / 2);
  },

  median: function median(mu/*, sigma*/) {
    return Math.exp(mu);
  },

  mode: function mode(mu, sigma) {
    return Math.exp(mu - sigma * sigma);
  },

  sample: function sample(mu, sigma) {
    return Math.exp(jStat.randn() * sigma + mu);
  },

  variance: function variance(mu, sigma) {
    return (Math.exp(sigma * sigma) - 1) * Math.exp(2 * mu + sigma * sigma);
  }
});



// extend noncentralt function with static methods
jStat.extend(jStat.noncentralt, {
  pdf: function pdf(x, dof, ncp) {
    var tol = 1e-14;
    if (Math.abs(ncp) < tol)  // ncp approx 0; use student-t
      return jStat.studentt.pdf(x, dof)

    if (Math.abs(x) < tol) {  // different formula for x == 0
      return Math.exp(jStat.gammaln((dof + 1) / 2) - ncp * ncp / 2 -
                      0.5 * Math.log(Math.PI * dof) - jStat.gammaln(dof / 2));
    }

    // formula for x != 0
    return dof / x *
        (jStat.noncentralt.cdf(x * Math.sqrt(1 + 2 / dof), dof+2, ncp) -
         jStat.noncentralt.cdf(x, dof, ncp));
  },

  cdf: function cdf(x, dof, ncp) {
    var tol = 1e-14;
    var min_iterations = 200;

    if (Math.abs(ncp) < tol)  // ncp approx 0; use student-t
      return jStat.studentt.cdf(x, dof);

    // turn negative x into positive and flip result afterwards
    var flip = false;
    if (x < 0) {
      flip = true;
      ncp = -ncp;
    }

    var prob = jStat.normal.cdf(-ncp, 0, 1);
    var value = tol + 1;
    // use value at last two steps to determine convergence
    var lastvalue = value;
    var y = x * x / (x * x + dof);
    var j = 0;
    var p = Math.exp(-ncp * ncp / 2);
    var q = Math.exp(-ncp * ncp / 2 - 0.5 * Math.log(2) -
                     jStat.gammaln(3 / 2)) * ncp;
    while (j < min_iterations || lastvalue > tol || value > tol) {
      lastvalue = value;
      if (j > 0) {
        p *= (ncp * ncp) / (2 * j);
        q *= (ncp * ncp) / (2 * (j + 1 / 2));
      }
      value = p * jStat.beta.cdf(y, j + 0.5, dof / 2) +
          q * jStat.beta.cdf(y, j+1, dof/2);
      prob += 0.5 * value;
      j++;
    }

    return flip ? (1 - prob) : prob;
  }
});


// extend normal function with static methods
jStat.extend(jStat.normal, {
  pdf: function pdf(x, mean, std) {
    return Math.exp(-0.5 * Math.log(2 * Math.PI) -
                    Math.log(std) - Math.pow(x - mean, 2) / (2 * std * std));
  },

  cdf: function cdf(x, mean, std) {
    return 0.5 * (1 + jStat.erf((x - mean) / Math.sqrt(2 * std * std)));
  },

  inv: function(p, mean, std) {
    return -1.41421356237309505 * std * jStat.erfcinv(2 * p) + mean;
  },

  mean : function(mean/*, std*/) {
    return mean;
  },

  median: function median(mean/*, std*/) {
    return mean;
  },

  mode: function (mean/*, std*/) {
    return mean;
  },

  sample: function sample(mean, std) {
    return jStat.randn() * std + mean;
  },

  variance : function(mean, std) {
    return std * std;
  }
});



// extend pareto function with static methods
jStat.extend(jStat.pareto, {
  pdf: function pdf(x, scale, shape) {
    if (x < scale)
      return 0;
    return (shape * Math.pow(scale, shape)) / Math.pow(x, shape + 1);
  },

  cdf: function cdf(x, scale, shape) {
    if (x < scale)
      return 0;
    return 1 - Math.pow(scale / x, shape);
  },

  inv: function inv(p, scale, shape) {
    return scale / Math.pow(1 - p, 1 / shape);
  },

  mean: function mean(scale, shape) {
    if (shape <= 1)
      return undefined;
    return (shape * Math.pow(scale, shape)) / (shape - 1);
  },

  median: function median(scale, shape) {
    return scale * (shape * Math.SQRT2);
  },

  mode: function mode(scale/*, shape*/) {
    return scale;
  },

  variance : function(scale, shape) {
    if (shape <= 2)
      return undefined;
    return (scale*scale * shape) / (Math.pow(shape - 1, 2) * (shape - 2));
  }
});



// extend studentt function with static methods
jStat.extend(jStat.studentt, {
  pdf: function pdf(x, dof) {
    dof = dof > 1e100 ? 1e100 : dof;
    return (1/(Math.sqrt(dof) * jStat.betafn(0.5, dof/2))) *
        Math.pow(1 + ((x * x) / dof), -((dof + 1) / 2));
  },

  cdf: function cdf(x, dof) {
    var dof2 = dof / 2;
    return jStat.ibeta((x + Math.sqrt(x * x + dof)) /
                       (2 * Math.sqrt(x * x + dof)), dof2, dof2);
  },

  inv: function(p, dof) {
    var x = jStat.ibetainv(2 * Math.min(p, 1 - p), 0.5 * dof, 0.5);
    x = Math.sqrt(dof * (1 - x) / x);
    return (p > 0.5) ? x : -x;
  },

  mean: function mean(dof) {
    return (dof > 1) ? 0 : undefined;
  },

  median: function median(/*dof*/) {
    return 0;
  },

  mode: function mode(/*dof*/) {
    return 0;
  },

  sample: function sample(dof) {
    return jStat.randn() * Math.sqrt(dof / (2 * jStat.randg(dof / 2)));
  },

  variance: function variance(dof) {
    return (dof  > 2) ? dof / (dof - 2) : (dof > 1) ? Infinity : undefined;
  }
});



// extend weibull function with static methods
jStat.extend(jStat.weibull, {
  pdf: function pdf(x, scale, shape) {
    if (x < 0 || scale < 0 || shape < 0)
      return 0;
    return (shape / scale) * Math.pow((x / scale), (shape - 1)) *
        Math.exp(-(Math.pow((x / scale), shape)));
  },

  cdf: function cdf(x, scale, shape) {
    return x < 0 ? 0 : 1 - Math.exp(-Math.pow((x / scale), shape));
  },

  inv: function(p, scale, shape) {
    return scale * Math.pow(-Math.log(1 - p), 1 / shape);
  },

  mean : function(scale, shape) {
    return scale * jStat.gammafn(1 + 1 / shape);
  },

  median: function median(scale, shape) {
    return scale * Math.pow(Math.log(2), 1 / shape);
  },

  mode: function mode(scale, shape) {
    if (shape <= 1)
      return 0;
    return scale * Math.pow((shape - 1) / shape, 1 / shape);
  },

  sample: function sample(scale, shape) {
    return scale * Math.pow(-Math.log(jStat._random_fn()), 1 / shape);
  },

  variance: function variance(scale, shape) {
    return scale * scale * jStat.gammafn(1 + 2 / shape) -
        Math.pow(jStat.weibull.mean(scale, shape), 2);
  }
});



// extend uniform function with static methods
jStat.extend(jStat.uniform, {
  pdf: function pdf(x, a, b) {
    return (x < a || x > b) ? 0 : 1 / (b - a);
  },

  cdf: function cdf(x, a, b) {
    if (x < a)
      return 0;
    else if (x < b)
      return (x - a) / (b - a);
    return 1;
  },

  inv: function(p, a, b) {
    return a + (p * (b - a));
  },

  mean: function mean(a, b) {
    return 0.5 * (a + b);
  },

  median: function median(a, b) {
    return jStat.mean(a, b);
  },

  mode: function mode(/*a, b*/) {
    throw new Error('mode is not yet implemented');
  },

  sample: function sample(a, b) {
    return (a / 2 + b / 2) + (b / 2 - a / 2) * (2 * jStat._random_fn() - 1);
  },

  variance: function variance(a, b) {
    return Math.pow(b - a, 2) / 12;
  }
});


// Got this from http://www.math.ucla.edu/~tom/distributions/binomial.html
function betinc(x, a, b, eps) {
  var a0 = 0;
  var b0 = 1;
  var a1 = 1;
  var b1 = 1;
  var m9 = 0;
  var a2 = 0;
  var c9;

  while (Math.abs((a1 - a2) / a1) > eps) {
    a2 = a1;
    c9 = -(a + m9) * (a + b + m9) * x / (a + 2 * m9) / (a + 2 * m9 + 1);
    a0 = a1 + c9 * a0;
    b0 = b1 + c9 * b0;
    m9 = m9 + 1;
    c9 = m9 * (b - m9) * x / (a + 2 * m9 - 1) / (a + 2 * m9);
    a1 = a0 + c9 * a1;
    b1 = b0 + c9 * b1;
    a0 = a0 / b1;
    b0 = b0 / b1;
    a1 = a1 / b1;
    b1 = 1;
  }

  return a1 / a;
}


// extend uniform function with static methods
jStat.extend(jStat.binomial, {
  pdf: function pdf(k, n, p) {
    return (p === 0 || p === 1) ?
      ((n * p) === k ? 1 : 0) :
      jStat.combination(n, k) * Math.pow(p, k) * Math.pow(1 - p, n - k);
  },

  cdf: function cdf(x, n, p) {
    var betacdf;
    var eps = 1e-10;

    if (x < 0)
      return 0;
    if (x >= n)
      return 1;
    if (p < 0 || p > 1 || n <= 0)
      return NaN;

    x = Math.floor(x);
    var z = p;
    var a = x + 1;
    var b = n - x;
    var s = a + b;
    var bt = Math.exp(jStat.gammaln(s) - jStat.gammaln(b) -
                      jStat.gammaln(a) + a * Math.log(z) + b * Math.log(1 - z));

    if (z < (a + 1) / (s + 2))
      betacdf = bt * betinc(z, a, b, eps);
    else
      betacdf = 1 - bt * betinc(1 - z, b, a, eps);

    return Math.round((1 - betacdf) * (1 / eps)) / (1 / eps);
  }
});



// extend uniform function with static methods
jStat.extend(jStat.negbin, {
  pdf: function pdf(k, r, p) {
    if (k !== k >>> 0)
      return false;
    if (k < 0)
      return 0;
    return jStat.combination(k + r - 1, r - 1) *
        Math.pow(1 - p, k) * Math.pow(p, r);
  },

  cdf: function cdf(x, r, p) {
    var sum = 0,
    k = 0;
    if (x < 0) return 0;
    for (; k <= x; k++) {
      sum += jStat.negbin.pdf(k, r, p);
    }
    return sum;
  }
});



// extend uniform function with static methods
jStat.extend(jStat.hypgeom, {
  pdf: function pdf(k, N, m, n) {
    // Hypergeometric PDF.

    // A simplification of the CDF algorithm below.

    // k = number of successes drawn
    // N = population size
    // m = number of successes in population
    // n = number of items drawn from population

    if(k !== k | 0) {
      return false;
    } else if(k < 0 || k < m - (N - n)) {
      // It's impossible to have this few successes drawn.
      return 0;
    } else if(k > n || k > m) {
      // It's impossible to have this many successes drawn.
      return 0;
    } else if (m * 2 > N) {
      // More than half the population is successes.

      if(n * 2 > N) {
        // More than half the population is sampled.

        return jStat.hypgeom.pdf(N - m - n + k, N, N - m, N - n)
      } else {
        // Half or less of the population is sampled.

        return jStat.hypgeom.pdf(n - k, N, N - m, n);
      }

    } else if(n * 2 > N) {
      // Half or less is successes.

      return jStat.hypgeom.pdf(m - k, N, m, N - n);

    } else if(m < n) {
      // We want to have the number of things sampled to be less than the
      // successes available. So swap the definitions of successful and sampled.
      return jStat.hypgeom.pdf(k, N, n, m);
    } else {
      // If we get here, half or less of the population was sampled, half or
      // less of it was successes, and we had fewer sampled things than
      // successes. Now we can do this complicated iterative algorithm in an
      // efficient way.

      // The basic premise of the algorithm is that we partially normalize our
      // intermediate product to keep it in a numerically good region, and then
      // finish the normalization at the end.

      // This variable holds the scaled probability of the current number of
      // successes.
      var scaledPDF = 1;

      // This keeps track of how much we have normalized.
      var samplesDone = 0;

      for(var i = 0; i < k; i++) {
        // For every possible number of successes up to that observed...

        while(scaledPDF > 1 && samplesDone < n) {
          // Intermediate result is growing too big. Apply some of the
          // normalization to shrink everything.

          scaledPDF *= 1 - (m / (N - samplesDone));

          // Say we've normalized by this sample already.
          samplesDone++;
        }

        // Work out the partially-normalized hypergeometric PDF for the next
        // number of successes
        scaledPDF *= (n - i) * (m - i) / ((i + 1) * (N - m - n + i + 1));
      }

      for(; samplesDone < n; samplesDone++) {
        // Apply all the rest of the normalization
        scaledPDF *= 1 - (m / (N - samplesDone));
      }

      // Bound answer sanely before returning.
      return Math.min(1, Math.max(0, scaledPDF));
    }
  },

  cdf: function cdf(x, N, m, n) {
    // Hypergeometric CDF.

    // This algorithm is due to Prof. Thomas S. Ferguson, <tom@math.ucla.edu>,
    // and comes from his hypergeometric test calculator at
    // <http://www.math.ucla.edu/~tom/distributions/Hypergeometric.html>.

    // x = number of successes drawn
    // N = population size
    // m = number of successes in population
    // n = number of items drawn from population

    if(x < 0 || x < m - (N - n)) {
      // It's impossible to have this few successes drawn or fewer.
      return 0;
    } else if(x >= n || x >= m) {
      // We will always have this many successes or fewer.
      return 1;
    } else if (m * 2 > N) {
      // More than half the population is successes.

      if(n * 2 > N) {
        // More than half the population is sampled.

        return jStat.hypgeom.cdf(N - m - n + x, N, N - m, N - n)
      } else {
        // Half or less of the population is sampled.

        return 1 - jStat.hypgeom.cdf(n - x - 1, N, N - m, n);
      }

    } else if(n * 2 > N) {
      // Half or less is successes.

      return 1 - jStat.hypgeom.cdf(m - x - 1, N, m, N - n);

    } else if(m < n) {
      // We want to have the number of things sampled to be less than the
      // successes available. So swap the definitions of successful and sampled.
      return jStat.hypgeom.cdf(x, N, n, m);
    } else {
      // If we get here, half or less of the population was sampled, half or
      // less of it was successes, and we had fewer sampled things than
      // successes. Now we can do this complicated iterative algorithm in an
      // efficient way.

      // The basic premise of the algorithm is that we partially normalize our
      // intermediate sum to keep it in a numerically good region, and then
      // finish the normalization at the end.

      // Holds the intermediate, scaled total CDF.
      var scaledCDF = 1;

      // This variable holds the scaled probability of the current number of
      // successes.
      var scaledPDF = 1;

      // This keeps track of how much we have normalized.
      var samplesDone = 0;

      for(var i = 0; i < x; i++) {
        // For every possible number of successes up to that observed...

        while(scaledCDF > 1 && samplesDone < n) {
          // Intermediate result is growing too big. Apply some of the
          // normalization to shrink everything.

          var factor = 1 - (m / (N - samplesDone));

          scaledPDF *= factor;
          scaledCDF *= factor;

          // Say we've normalized by this sample already.
          samplesDone++;
        }

        // Work out the partially-normalized hypergeometric PDF for the next
        // number of successes
        scaledPDF *= (n - i) * (m - i) / ((i + 1) * (N - m - n + i + 1));

        // Add to the CDF answer.
        scaledCDF += scaledPDF;
      }

      for(; samplesDone < n; samplesDone++) {
        // Apply all the rest of the normalization
        scaledCDF *= 1 - (m / (N - samplesDone));
      }

      // Bound answer sanely before returning.
      return Math.min(1, Math.max(0, scaledCDF));
    }
  }
});



// extend uniform function with static methods
jStat.extend(jStat.poisson, {
  pdf: function pdf(k, l) {
    if (l < 0 || (k % 1) !== 0 || k < 0) {
      return 0;
    }

    return Math.pow(l, k) * Math.exp(-l) / jStat.factorial(k);
  },

  cdf: function cdf(x, l) {
    var sumarr = [],
    k = 0;
    if (x < 0) return 0;
    for (; k <= x; k++) {
      sumarr.push(jStat.poisson.pdf(k, l));
    }
    return jStat.sum(sumarr);
  },

  mean : function(l) {
    return l;
  },

  variance : function(l) {
    return l;
  },

  sampleSmall: function sampleSmall(l) {
    var p = 1, k = 0, L = Math.exp(-l);
    do {
      k++;
      p *= jStat._random_fn();
    } while (p > L);
    return k - 1;
  },

  sampleLarge: function sampleLarge(l) {
    var lam = l;
    var k;
    var U, V, slam, loglam, a, b, invalpha, vr, us;

    slam = Math.sqrt(lam);
    loglam = Math.log(lam);
    b = 0.931 + 2.53 * slam;
    a = -0.059 + 0.02483 * b;
    invalpha = 1.1239 + 1.1328 / (b - 3.4);
    vr = 0.9277 - 3.6224 / (b - 2);

    while (1) {
      U = Math.random() - 0.5;
      V = Math.random();
      us = 0.5 - Math.abs(U);
      k = Math.floor((2 * a / us + b) * U + lam + 0.43);
      if ((us >= 0.07) && (V <= vr)) {
          return k;
      }
      if ((k < 0) || ((us < 0.013) && (V > us))) {
          continue;
      }
      /* log(V) == log(0.0) ok here */
      /* if U==0.0 so that us==0.0, log is ok since always returns */
      if ((Math.log(V) + Math.log(invalpha) - Math.log(a / (us * us) + b)) <= (-lam + k * loglam - jStat.loggam(k + 1))) {
          return k;
      }
    }
  },

  sample: function sample(l) {
    if (l < 10)
      return this.sampleSmall(l);
    else
      return this.sampleLarge(l);
  }
});

// extend triangular function with static methods
jStat.extend(jStat.triangular, {
  pdf: function pdf(x, a, b, c) {
    if (b <= a || c < a || c > b) {
      return NaN;
    } else {
      if (x < a || x > b) {
        return 0;
      } else if (x < c) {
          return (2 * (x - a)) / ((b - a) * (c - a));
      } else if (x === c) {
          return (2 / (b - a));
      } else { // x > c
          return (2 * (b - x)) / ((b - a) * (b - c));
      }
    }
  },

  cdf: function cdf(x, a, b, c) {
    if (b <= a || c < a || c > b)
      return NaN;
    if (x <= a)
      return 0;
    else if (x >= b)
      return 1;
    if (x <= c)
      return Math.pow(x - a, 2) / ((b - a) * (c - a));
    else // x > c
      return 1 - Math.pow(b - x, 2) / ((b - a) * (b - c));
  },

  inv: function inv(p, a, b, c) {
    if (b <= a || c < a || c > b) {
      return NaN;
    } else {
      if (p <= ((c - a) / (b - a))) {
        return a + (b - a) * Math.sqrt(p * ((c - a) / (b - a)));
      } else { // p > ((c - a) / (b - a))
        return a + (b - a) * (1 - Math.sqrt((1 - p) * (1 - ((c - a) / (b - a)))));
      }
    }
  },

  mean: function mean(a, b, c) {
    return (a + b + c) / 3;
  },

  median: function median(a, b, c) {
    if (c <= (a + b) / 2) {
      return b - Math.sqrt((b - a) * (b - c)) / Math.sqrt(2);
    } else if (c > (a + b) / 2) {
      return a + Math.sqrt((b - a) * (c - a)) / Math.sqrt(2);
    }
  },

  mode: function mode(a, b, c) {
    return c;
  },

  sample: function sample(a, b, c) {
    var u = jStat._random_fn();
    if (u < ((c - a) / (b - a)))
      return a + Math.sqrt(u * (b - a) * (c - a))
    return b - Math.sqrt((1 - u) * (b - a) * (b - c));
  },

  variance: function variance(a, b, c) {
    return (a * a + b * b + c * c - a * b - a * c - b * c) / 18;
  }
});


// extend arcsine function with static methods
jStat.extend(jStat.arcsine, {
  pdf: function pdf(x, a, b) {
    if (b <= a) return NaN;

    return (x <= a || x >= b) ? 0 :
      (2 / Math.PI) *
        Math.pow(Math.pow(b - a, 2) -
                  Math.pow(2 * x - a - b, 2), -0.5);
  },

  cdf: function cdf(x, a, b) {
    if (x < a)
      return 0;
    else if (x < b)
      return (2 / Math.PI) * Math.asin(Math.sqrt((x - a)/(b - a)));
    return 1;
  },

  inv: function(p, a, b) {
    return a + (0.5 - 0.5 * Math.cos(Math.PI * p)) * (b - a);
  },

  mean: function mean(a, b) {
    if (b <= a) return NaN;
    return (a + b) / 2;
  },

  median: function median(a, b) {
    if (b <= a) return NaN;
    return (a + b) / 2;
  },

  mode: function mode(/*a, b*/) {
    throw new Error('mode is not yet implemented');
  },

  sample: function sample(a, b) {
    return ((a + b) / 2) + ((b - a) / 2) *
      Math.sin(2 * Math.PI * jStat.uniform.sample(0, 1));
  },

  variance: function variance(a, b) {
    if (b <= a) return NaN;
    return Math.pow(b - a, 2) / 8;
  }
});


function laplaceSign(x) { return x / Math.abs(x); }

jStat.extend(jStat.laplace, {
  pdf: function pdf(x, mu, b) {
    return (b <= 0) ? 0 : (Math.exp(-Math.abs(x - mu) / b)) / (2 * b);
  },

  cdf: function cdf(x, mu, b) {
    if (b <= 0) { return 0; }

    if(x < mu) {
      return 0.5 * Math.exp((x - mu) / b);
    } else {
      return 1 - 0.5 * Math.exp(- (x - mu) / b);
    }
  },

  mean: function(mu/*, b*/) {
    return mu;
  },

  median: function(mu/*, b*/) {
    return mu;
  },

  mode: function(mu/*, b*/) {
    return mu;
  },

  variance: function(mu, b) {
    return 2 * b * b;
  },

  sample: function sample(mu, b) {
    var u = jStat._random_fn() - 0.5;

    return mu - (b * laplaceSign(u) * Math.log(1 - (2 * Math.abs(u))));
  }
});

function tukeyWprob(w, rr, cc) {
  var nleg = 12;
  var ihalf = 6;

  var C1 = -30;
  var C2 = -50;
  var C3 = 60;
  var bb   = 8;
  var wlar = 3;
  var wincr1 = 2;
  var wincr2 = 3;
  var xleg = [
    0.981560634246719250690549090149,
    0.904117256370474856678465866119,
    0.769902674194304687036893833213,
    0.587317954286617447296702418941,
    0.367831498998180193752691536644,
    0.125233408511468915472441369464
  ];
  var aleg = [
    0.047175336386511827194615961485,
    0.106939325995318430960254718194,
    0.160078328543346226334652529543,
    0.203167426723065921749064455810,
    0.233492536538354808760849898925,
    0.249147045813402785000562436043
  ];

  var qsqz = w * 0.5;

  // if w >= 16 then the integral lower bound (occurs for c=20)
  // is 0.99999999999995 so return a value of 1.

  if (qsqz >= bb)
    return 1.0;

  // find (f(w/2) - 1) ^ cc
  // (first term in integral of hartley's form).

  var pr_w = 2 * jStat.normal.cdf(qsqz, 0, 1, 1, 0) - 1; // erf(qsqz / M_SQRT2)
  // if pr_w ^ cc < 2e-22 then set pr_w = 0
  if (pr_w >= Math.exp(C2 / cc))
    pr_w = Math.pow(pr_w, cc);
  else
    pr_w = 0.0;

  // if w is large then the second component of the
  // integral is small, so fewer intervals are needed.

  var wincr;
  if (w > wlar)
    wincr = wincr1;
  else
    wincr = wincr2;

  // find the integral of second term of hartley's form
  // for the integral of the range for equal-length
  // intervals using legendre quadrature.  limits of
  // integration are from (w/2, 8).  two or three
  // equal-length intervals are used.

  // blb and bub are lower and upper limits of integration.

  var blb = qsqz;
  var binc = (bb - qsqz) / wincr;
  var bub = blb + binc;
  var einsum = 0.0;

  // integrate over each interval

  var cc1 = cc - 1.0;
  for (var wi = 1; wi <= wincr; wi++) {
    var elsum = 0.0;
    var a = 0.5 * (bub + blb);

    // legendre quadrature with order = nleg

    var b = 0.5 * (bub - blb);

    for (var jj = 1; jj <= nleg; jj++) {
      var j, xx;
      if (ihalf < jj) {
        j = (nleg - jj) + 1;
        xx = xleg[j-1];
      } else {
        j = jj;
        xx = -xleg[j-1];
      }
      var c = b * xx;
      var ac = a + c;

      // if exp(-qexpo/2) < 9e-14,
      // then doesn't contribute to integral

      var qexpo = ac * ac;
      if (qexpo > C3)
        break;

      var pplus = 2 * jStat.normal.cdf(ac, 0, 1, 1, 0);
      var pminus= 2 * jStat.normal.cdf(ac, w, 1, 1, 0);

      // if rinsum ^ (cc-1) < 9e-14,
      // then doesn't contribute to integral

      var rinsum = (pplus * 0.5) - (pminus * 0.5);
      if (rinsum >= Math.exp(C1 / cc1)) {
        rinsum = (aleg[j-1] * Math.exp(-(0.5 * qexpo))) * Math.pow(rinsum, cc1);
        elsum += rinsum;
      }
    }
    elsum *= (((2.0 * b) * cc) / Math.sqrt(2 * Math.PI));
    einsum += elsum;
    blb = bub;
    bub += binc;
  }

  // if pr_w ^ rr < 9e-14, then return 0
  pr_w += einsum;
  if (pr_w <= Math.exp(C1 / rr))
    return 0;

  pr_w = Math.pow(pr_w, rr);
  if (pr_w >= 1) // 1 was iMax was eps
    return 1;
  return pr_w;
}

function tukeyQinv(p, c, v) {
  var p0 = 0.322232421088;
  var q0 = 0.993484626060e-01;
  var p1 = -1.0;
  var q1 = 0.588581570495;
  var p2 = -0.342242088547;
  var q2 = 0.531103462366;
  var p3 = -0.204231210125;
  var q3 = 0.103537752850;
  var p4 = -0.453642210148e-04;
  var q4 = 0.38560700634e-02;
  var c1 = 0.8832;
  var c2 = 0.2368;
  var c3 = 1.214;
  var c4 = 1.208;
  var c5 = 1.4142;
  var vmax = 120.0;

  var ps = 0.5 - 0.5 * p;
  var yi = Math.sqrt(Math.log(1.0 / (ps * ps)));
  var t = yi + (((( yi * p4 + p3) * yi + p2) * yi + p1) * yi + p0)
     / (((( yi * q4 + q3) * yi + q2) * yi + q1) * yi + q0);
  if (v < vmax) t += (t * t * t + t) / v / 4.0;
  var q = c1 - c2 * t;
  if (v < vmax) q += -c3 / v + c4 * t / v;
  return t * (q * Math.log(c - 1.0) + c5);
}

jStat.extend(jStat.tukey, {
  cdf: function cdf(q, nmeans, df) {
    // Identical implementation as the R ptukey() function as of commit 68947
    var rr = 1;
    var cc = nmeans;

    var nlegq = 16;
    var ihalfq = 8;

    var eps1 = -30.0;
    var eps2 = 1.0e-14;
    var dhaf  = 100.0;
    var dquar = 800.0;
    var deigh = 5000.0;
    var dlarg = 25000.0;
    var ulen1 = 1.0;
    var ulen2 = 0.5;
    var ulen3 = 0.25;
    var ulen4 = 0.125;
    var xlegq = [
      0.989400934991649932596154173450,
      0.944575023073232576077988415535,
      0.865631202387831743880467897712,
      0.755404408355003033895101194847,
      0.617876244402643748446671764049,
      0.458016777657227386342419442984,
      0.281603550779258913230460501460,
      0.950125098376374401853193354250e-1
    ];
    var alegq = [
      0.271524594117540948517805724560e-1,
      0.622535239386478928628438369944e-1,
      0.951585116824927848099251076022e-1,
      0.124628971255533872052476282192,
      0.149595988816576732081501730547,
      0.169156519395002538189312079030,
      0.182603415044923588866763667969,
      0.189450610455068496285396723208
    ];

    if (q <= 0)
      return 0;

    // df must be > 1
    // there must be at least two values

    if (df < 2 || rr < 1 || cc < 2) return NaN;

    if (!Number.isFinite(q))
      return 1;

    if (df > dlarg)
      return tukeyWprob(q, rr, cc);

    // calculate leading constant

    var f2 = df * 0.5;
    var f2lf = ((f2 * Math.log(df)) - (df * Math.log(2))) - jStat.gammaln(f2);
    var f21 = f2 - 1.0;

    // integral is divided into unit, half-unit, quarter-unit, or
    // eighth-unit length intervals depending on the value of the
    // degrees of freedom.

    var ff4 = df * 0.25;
    var ulen;
    if      (df <= dhaf)  ulen = ulen1;
    else if (df <= dquar) ulen = ulen2;
    else if (df <= deigh) ulen = ulen3;
    else                  ulen = ulen4;

    f2lf += Math.log(ulen);

    // integrate over each subinterval

    var ans = 0.0;

    for (var i = 1; i <= 50; i++) {
      var otsum = 0.0;

      // legendre quadrature with order = nlegq
      // nodes (stored in xlegq) are symmetric around zero.

      var twa1 = (2 * i - 1) * ulen;

      for (var jj = 1; jj <= nlegq; jj++) {
        var j, t1;
        if (ihalfq < jj) {
          j = jj - ihalfq - 1;
          t1 = (f2lf + (f21 * Math.log(twa1 + (xlegq[j] * ulen))))
              - (((xlegq[j] * ulen) + twa1) * ff4);
        } else {
          j = jj - 1;
          t1 = (f2lf + (f21 * Math.log(twa1 - (xlegq[j] * ulen))))
              + (((xlegq[j] * ulen) - twa1) * ff4);
        }

        // if exp(t1) < 9e-14, then doesn't contribute to integral
        var qsqz;
        if (t1 >= eps1) {
          if (ihalfq < jj) {
            qsqz = q * Math.sqrt(((xlegq[j] * ulen) + twa1) * 0.5);
          } else {
            qsqz = q * Math.sqrt(((-(xlegq[j] * ulen)) + twa1) * 0.5);
          }

          // call wprob to find integral of range portion

          var wprb = tukeyWprob(qsqz, rr, cc);
          var rotsum = (wprb * alegq[j]) * Math.exp(t1);
          otsum += rotsum;
        }
        // end legendre integral for interval i
        // L200:
      }

      // if integral for interval i < 1e-14, then stop.
      // However, in order to avoid small area under left tail,
      // at least  1 / ulen  intervals are calculated.
      if (i * ulen >= 1.0 && otsum <= eps2)
        break;

      // end of interval i
      // L330:

      ans += otsum;
    }

    if (otsum > eps2) { // not converged
      throw new Error('tukey.cdf failed to converge');
    }
    if (ans > 1)
      ans = 1;
    return ans;
  },

  inv: function(p, nmeans, df) {
    // Identical implementation as the R qtukey() function as of commit 68947
    var rr = 1;
    var cc = nmeans;

    var eps = 0.0001;
    var maxiter = 50;

    // df must be > 1 ; there must be at least two values
    if (df < 2 || rr < 1 || cc < 2) return NaN;

    if (p < 0 || p > 1) return NaN;
    if (p === 0) return 0;
    if (p === 1) return Infinity;

    // Initial value

    var x0 = tukeyQinv(p, cc, df);

    // Find prob(value < x0)

    var valx0 = jStat.tukey.cdf(x0, nmeans, df) - p;

    // Find the second iterate and prob(value < x1).
    // If the first iterate has probability value
    // exceeding p then second iterate is 1 less than
    // first iterate; otherwise it is 1 greater.

    var x1;
    if (valx0 > 0.0)
      x1 = Math.max(0.0, x0 - 1.0);
    else
      x1 = x0 + 1.0;
    var valx1 = jStat.tukey.cdf(x1, nmeans, df) - p;

    // Find new iterate

    var ans;
    for(var iter = 1; iter < maxiter; iter++) {
      ans = x1 - ((valx1 * (x1 - x0)) / (valx1 - valx0));
      valx0 = valx1;

      // New iterate must be >= 0

      x0 = x1;
      if (ans < 0.0) {
        ans = 0.0;
        valx1 = -p;
      }
      // Find prob(value < new iterate)

      valx1 = jStat.tukey.cdf(ans, nmeans, df) - p;
      x1 = ans;

      // If the difference between two successive
      // iterates is less than eps, stop

      var xabs = Math.abs(x1 - x0);
      if (xabs < eps)
        return ans;
    }

    throw new Error('tukey.inv failed to converge');
  }
});

}(jStat, Math));
/* Provides functions for the solution of linear system of equations, integration, extrapolation,
 * interpolation, eigenvalue problems, differential equations and PCA analysis. */

(function(jStat, Math) {

var push = Array.prototype.push;
var isArray = jStat.utils.isArray;

function isUsable(arg) {
  return isArray(arg) || arg instanceof jStat;
}

jStat.extend({

  // add a vector/matrix to a vector/matrix or scalar
  add: function add(arr, arg) {
    // check if arg is a vector or scalar
    if (isUsable(arg)) {
      if (!isUsable(arg[0])) arg = [ arg ];
      return jStat.map(arr, function(value, row, col) {
        return value + arg[row][col];
      });
    }
    return jStat.map(arr, function(value) { return value + arg; });
  },

  // subtract a vector or scalar from the vector
  subtract: function subtract(arr, arg) {
    // check if arg is a vector or scalar
    if (isUsable(arg)) {
      if (!isUsable(arg[0])) arg = [ arg ];
      return jStat.map(arr, function(value, row, col) {
        return value - arg[row][col] || 0;
      });
    }
    return jStat.map(arr, function(value) { return value - arg; });
  },

  // matrix division
  divide: function divide(arr, arg) {
    if (isUsable(arg)) {
      if (!isUsable(arg[0])) arg = [ arg ];
      return jStat.multiply(arr, jStat.inv(arg));
    }
    return jStat.map(arr, function(value) { return value / arg; });
  },

  // matrix multiplication
  multiply: function multiply(arr, arg) {
    var row, col, nrescols, sum, nrow, ncol, res, rescols;
    // eg: arr = 2 arg = 3 -> 6 for res[0][0] statement closure
    if (arr.length === undefined && arg.length === undefined) {
      return arr * arg;
    }
    nrow = arr.length,
    ncol = arr[0].length,
    res = jStat.zeros(nrow, nrescols = (isUsable(arg)) ? arg[0].length : ncol),
    rescols = 0;
    if (isUsable(arg)) {
      for (; rescols < nrescols; rescols++) {
        for (row = 0; row < nrow; row++) {
          sum = 0;
          for (col = 0; col < ncol; col++)
          sum += arr[row][col] * arg[col][rescols];
          res[row][rescols] = sum;
        }
      }
      return (nrow === 1 && rescols === 1) ? res[0][0] : res;
    }
    return jStat.map(arr, function(value) { return value * arg; });
  },

  // outer([1,2,3],[4,5,6])
  // ===
  // [[1],[2],[3]] times [[4,5,6]]
  // ->
  // [[4,5,6],[8,10,12],[12,15,18]]
  outer:function outer(A, B) {
    return jStat.multiply(A.map(function(t){ return [t] }), [B]);
  },


  // Returns the dot product of two matricies
  dot: function dot(arr, arg) {
    if (!isUsable(arr[0])) arr = [ arr ];
    if (!isUsable(arg[0])) arg = [ arg ];
    // convert column to row vector
    var left = (arr[0].length === 1 && arr.length !== 1) ? jStat.transpose(arr) : arr,
    right = (arg[0].length === 1 && arg.length !== 1) ? jStat.transpose(arg) : arg,
    res = [],
    row = 0,
    nrow = left.length,
    ncol = left[0].length,
    sum, col;
    for (; row < nrow; row++) {
      res[row] = [];
      sum = 0;
      for (col = 0; col < ncol; col++)
      sum += left[row][col] * right[row][col];
      res[row] = sum;
    }
    return (res.length === 1) ? res[0] : res;
  },

  // raise every element by a scalar
  pow: function pow(arr, arg) {
    return jStat.map(arr, function(value) { return Math.pow(value, arg); });
  },

  // exponentiate every element
  exp: function exp(arr) {
    return jStat.map(arr, function(value) { return Math.exp(value); });
  },

  // generate the natural log of every element
  log: function exp(arr) {
    return jStat.map(arr, function(value) { return Math.log(value); });
  },

  // generate the absolute values of the vector
  abs: function abs(arr) {
    return jStat.map(arr, function(value) { return Math.abs(value); });
  },

  // computes the p-norm of the vector
  // In the case that a matrix is passed, uses the first row as the vector
  norm: function norm(arr, p) {
    var nnorm = 0,
    i = 0;
    // check the p-value of the norm, and set for most common case
    if (isNaN(p)) p = 2;
    // check if multi-dimensional array, and make vector correction
    if (isUsable(arr[0])) arr = arr[0];
    // vector norm
    for (; i < arr.length; i++) {
      nnorm += Math.pow(Math.abs(arr[i]), p);
    }
    return Math.pow(nnorm, 1 / p);
  },

  // computes the angle between two vectors in rads
  // In case a matrix is passed, this uses the first row as the vector
  angle: function angle(arr, arg) {
    return Math.acos(jStat.dot(arr, arg) / (jStat.norm(arr) * jStat.norm(arg)));
  },

  // augment one matrix by another
  // Note: this function returns a matrix, not a jStat object
  aug: function aug(a, b) {
    var newarr = [];
    var i;
    for (i = 0; i < a.length; i++) {
      newarr.push(a[i].slice());
    }
    for (i = 0; i < newarr.length; i++) {
      push.apply(newarr[i], b[i]);
    }
    return newarr;
  },

  // The inv() function calculates the inverse of a matrix
  // Create the inverse by augmenting the matrix by the identity matrix of the
  // appropriate size, and then use G-J elimination on the augmented matrix.
  inv: function inv(a) {
    var rows = a.length;
    var cols = a[0].length;
    var b = jStat.identity(rows, cols);
    var c = jStat.gauss_jordan(a, b);
    var result = [];
    var i = 0;
    var j;

    //We need to copy the inverse portion to a new matrix to rid G-J artifacts
    for (; i < rows; i++) {
      result[i] = [];
      for (j = cols; j < c[0].length; j++)
        result[i][j - cols] = c[i][j];
    }
    return result;
  },

  // calculate the determinant of a matrix
  det: function det(a) {
    if (a.length === 2) {
      return a[0][0] * a[1][1] - a[0][1] * a[1][0];
    }

    var determinant = 0;
    for (var i = 0; i < a.length; i++) {
      // build a sub matrix without column `i`
      var submatrix = [];
      for (var row = 1; row < a.length; row++) {
        submatrix[row - 1] = [];
        for (var col = 0; col < a.length; col++) {
          if (col < i) {
            submatrix[row - 1][col] = a[row][col];
          } else if (col > i) {
            submatrix[row - 1][col - 1] = a[row][col];
          }
        }
      }

      // alternate between + and - between determinants
      var sign = i % 2 ? -1 : 1;
      determinant += det(submatrix) * a[0][i] * sign;
    }

    return determinant
  },

  gauss_elimination: function gauss_elimination(a, b) {
    var i = 0,
    j = 0,
    n = a.length,
    m = a[0].length,
    factor = 1,
    sum = 0,
    x = [],
    maug, pivot, temp, k;
    a = jStat.aug(a, b);
    maug = a[0].length;
    for(i = 0; i < n; i++) {
      pivot = a[i][i];
      j = i;
      for (k = i + 1; k < m; k++) {
        if (pivot < Math.abs(a[k][i])) {
          pivot = a[k][i];
          j = k;
        }
      }
      if (j != i) {
        for(k = 0; k < maug; k++) {
          temp = a[i][k];
          a[i][k] = a[j][k];
          a[j][k] = temp;
        }
      }
      for (j = i + 1; j < n; j++) {
        factor = a[j][i] / a[i][i];
        for(k = i; k < maug; k++) {
          a[j][k] = a[j][k] - factor * a[i][k];
        }
      }
    }
    for (i = n - 1; i >= 0; i--) {
      sum = 0;
      for (j = i + 1; j<= n - 1; j++) {
        sum = sum + x[j] * a[i][j];
      }
      x[i] =(a[i][maug - 1] - sum) / a[i][i];
    }
    return x;
  },

  gauss_jordan: function gauss_jordan(a, b) {
    var m = jStat.aug(a, b);
    var h = m.length;
    var w = m[0].length;
    var c = 0;
    var x, y, y2;
    // find max pivot
    for (y = 0; y < h; y++) {
      var maxrow = y;
      for (y2 = y+1; y2 < h; y2++) {
        if (Math.abs(m[y2][y]) > Math.abs(m[maxrow][y]))
          maxrow = y2;
      }
      var tmp = m[y];
      m[y] = m[maxrow];
      m[maxrow] = tmp
      for (y2 = y+1; y2 < h; y2++) {
        c = m[y2][y] / m[y][y];
        for (x = y; x < w; x++) {
          m[y2][x] -= m[y][x] * c;
        }
      }
    }
    // backsubstitute
    for (y = h-1; y >= 0; y--) {
      c = m[y][y];
      for (y2 = 0; y2 < y; y2++) {
        for (x = w-1; x > y-1; x--) {
          m[y2][x] -= m[y][x] * m[y2][y] / c;
        }
      }
      m[y][y] /= c;
      for (x = h; x < w; x++) {
        m[y][x] /= c;
      }
    }
    return m;
  },

  // solve equation
  // Ax=b
  // A is upper triangular matrix
  // A=[[1,2,3],[0,4,5],[0,6,7]]
  // b=[1,2,3]
  // triaUpSolve(A,b) // -> [2.666,0.1666,1.666]
  // if you use matrix style
  // A=[[1,2,3],[0,4,5],[0,6,7]]
  // b=[[1],[2],[3]]
  // will return [[2.666],[0.1666],[1.666]]
  triaUpSolve: function triaUpSolve(A, b) {
    var size = A[0].length;
    var x = jStat.zeros(1, size)[0];
    var parts;
    var matrix_mode = false;

    if (b[0].length != undefined) {
      b = b.map(function(i){ return i[0] });
      matrix_mode = true;
    }

    jStat.arange(size - 1, -1, -1).forEach(function(i) {
      parts = jStat.arange(i + 1, size).map(function(j) {
        return x[j] * A[i][j];
      });
      x[i] = (b[i] - jStat.sum(parts)) / A[i][i];
    });

    if (matrix_mode)
      return x.map(function(i){ return [i] });
    return x;
  },

  triaLowSolve: function triaLowSolve(A, b) {
    // like to triaUpSolve but A is lower triangular matrix
    var size = A[0].length;
    var x = jStat.zeros(1, size)[0];
    var parts;

    var matrix_mode=false;
    if (b[0].length != undefined) {
      b = b.map(function(i){ return i[0] });
      matrix_mode = true;
    }

    jStat.arange(size).forEach(function(i) {
      parts = jStat.arange(i).map(function(j) {
        return A[i][j] * x[j];
      });
      x[i] = (b[i] - jStat.sum(parts)) / A[i][i];
    })

    if (matrix_mode)
      return x.map(function(i){ return [i] });
    return x;
  },


  // A -> [L,U]
  // A=LU
  // L is lower triangular matrix
  // U is upper triangular matrix
  lu: function lu(A) {
    var size = A.length;
    //var L=jStat.diagonal(jStat.ones(1,size)[0]);
    var L = jStat.identity(size);
    var R = jStat.zeros(A.length, A[0].length);
    var parts;
    jStat.arange(size).forEach(function(t) {
      R[0][t] = A[0][t];
    });
    jStat.arange(1, size).forEach(function(l) {
      jStat.arange(l).forEach(function(i) {
        parts = jStat.arange(i).map(function(jj) {
          return L[l][jj] * R[jj][i];
        });
        L[l][i] = (A[l][i] - jStat.sum(parts)) / R[i][i];
      });
      jStat.arange(l, size).forEach(function(j) {
        parts = jStat.arange(l).map(function(jj) {
          return L[l][jj] * R[jj][j];
        });
        R[l][j] = A[parts.length][j] - jStat.sum(parts);
      });
    });
    return [L, R];
  },

  // A -> T
  // A=TT'
  // T is lower triangular matrix
  cholesky: function cholesky(A) {
    var size = A.length;
    var T = jStat.zeros(A.length, A[0].length);
    var parts;
    jStat.arange(size).forEach(function(i) {
      parts = jStat.arange(i).map(function(t) {
        return Math.pow(T[i][t],2);
      });
      T[i][i] = Math.sqrt(A[i][i] - jStat.sum(parts));
      jStat.arange(i + 1, size).forEach(function(j) {
        parts = jStat.arange(i).map(function(t) {
          return T[i][t] * T[j][t];
        });
        T[j][i] = (A[i][j] - jStat.sum(parts)) / T[i][i];
      });
    });
    return T;
  },


  gauss_jacobi: function gauss_jacobi(a, b, x, r) {
    var i = 0;
    var j = 0;
    var n = a.length;
    var l = [];
    var u = [];
    var d = [];
    var xv, c, h, xk;
    for (; i < n; i++) {
      l[i] = [];
      u[i] = [];
      d[i] = [];
      for (j = 0; j < n; j++) {
        if (i > j) {
          l[i][j] = a[i][j];
          u[i][j] = d[i][j] = 0;
        } else if (i < j) {
          u[i][j] = a[i][j];
          l[i][j] = d[i][j] = 0;
        } else {
          d[i][j] = a[i][j];
          l[i][j] = u[i][j] = 0;
        }
      }
    }
    h = jStat.multiply(jStat.multiply(jStat.inv(d), jStat.add(l, u)), -1);
    c = jStat.multiply(jStat.inv(d), b);
    xv = x;
    xk = jStat.add(jStat.multiply(h, x), c);
    i = 2;
    while (Math.abs(jStat.norm(jStat.subtract(xk,xv))) > r) {
      xv = xk;
      xk = jStat.add(jStat.multiply(h, xv), c);
      i++;
    }
    return xk;
  },

  gauss_seidel: function gauss_seidel(a, b, x, r) {
    var i = 0;
    var n = a.length;
    var l = [];
    var u = [];
    var d = [];
    var j, xv, c, h, xk;
    for (; i < n; i++) {
      l[i] = [];
      u[i] = [];
      d[i] = [];
      for (j = 0; j < n; j++) {
        if (i > j) {
          l[i][j] = a[i][j];
          u[i][j] = d[i][j] = 0;
        } else if (i < j) {
          u[i][j] = a[i][j];
          l[i][j] = d[i][j] = 0;
        } else {
          d[i][j] = a[i][j];
          l[i][j] = u[i][j] = 0;
        }
      }
    }
    h = jStat.multiply(jStat.multiply(jStat.inv(jStat.add(d, l)), u), -1);
    c = jStat.multiply(jStat.inv(jStat.add(d, l)), b);
    xv = x;
    xk = jStat.add(jStat.multiply(h, x), c);
    i = 2;
    while (Math.abs(jStat.norm(jStat.subtract(xk, xv))) > r) {
      xv = xk;
      xk = jStat.add(jStat.multiply(h, xv), c);
      i = i + 1;
    }
    return xk;
  },

  SOR: function SOR(a, b, x, r, w) {
    var i = 0;
    var n = a.length;
    var l = [];
    var u = [];
    var d = [];
    var j, xv, c, h, xk;
    for (; i < n; i++) {
      l[i] = [];
      u[i] = [];
      d[i] = [];
      for (j = 0; j < n; j++) {
        if (i > j) {
          l[i][j] = a[i][j];
          u[i][j] = d[i][j] = 0;
        } else if (i < j) {
          u[i][j] = a[i][j];
          l[i][j] = d[i][j] = 0;
        } else {
          d[i][j] = a[i][j];
          l[i][j] = u[i][j] = 0;
        }
      }
    }
    h = jStat.multiply(jStat.inv(jStat.add(d, jStat.multiply(l, w))),
                       jStat.subtract(jStat.multiply(d, 1 - w),
                                      jStat.multiply(u, w)));
    c = jStat.multiply(jStat.multiply(jStat.inv(jStat.add(d,
        jStat.multiply(l, w))), b), w);
    xv = x;
    xk = jStat.add(jStat.multiply(h, x), c);
    i = 2;
    while (Math.abs(jStat.norm(jStat.subtract(xk, xv))) > r) {
      xv = xk;
      xk = jStat.add(jStat.multiply(h, xv), c);
      i++;
    }
    return xk;
  },

  householder: function householder(a) {
    var m = a.length;
    var n = a[0].length;
    var i = 0;
    var w = [];
    var p = [];
    var alpha, r, k, j, factor;
    for (; i < m - 1; i++) {
      alpha = 0;
      for (j = i + 1; j < n; j++)
      alpha += (a[j][i] * a[j][i]);
      factor = (a[i + 1][i] > 0) ? -1 : 1;
      alpha = factor * Math.sqrt(alpha);
      r = Math.sqrt((((alpha * alpha) - a[i + 1][i] * alpha) / 2));
      w = jStat.zeros(m, 1);
      w[i + 1][0] = (a[i + 1][i] - alpha) / (2 * r);
      for (k = i + 2; k < m; k++) w[k][0] = a[k][i] / (2 * r);
      p = jStat.subtract(jStat.identity(m, n),
          jStat.multiply(jStat.multiply(w, jStat.transpose(w)), 2));
      a = jStat.multiply(p, jStat.multiply(a, p));
    }
    return a;
  },

  // A -> [Q,R]
  // Q is orthogonal matrix
  // R is upper triangular
  QR: (function() {
    // x -> Q
    // find a orthogonal matrix Q st.
    // Qx=y
    // y is [||x||,0,0,...]

    // quick ref
    var sum   = jStat.sum;
    var range = jStat.arange;

    function qr2(x) {
      // quick impletation
      // https://www.stat.wisc.edu/~larget/math496/qr.html

      var n = x.length;
      var p = x[0].length;

      var r = jStat.zeros(p, p);
      x = jStat.copy(x);

      var i,j,k;
      for(j = 0; j < p; j++){
        r[j][j] = Math.sqrt(sum(range(n).map(function(i){
          return x[i][j] * x[i][j];
        })));
        for(i = 0; i < n; i++){
          x[i][j] = x[i][j] / r[j][j];
        }
        for(k = j+1; k < p; k++){
          r[j][k] = sum(range(n).map(function(i){
            return x[i][j] * x[i][k];
          }));
          for(i = 0; i < n; i++){
            x[i][k] = x[i][k] - x[i][j]*r[j][k];
          }
        }
      }
      return [x, r];
    }

    return qr2;
  }()),

  lstsq: (function() {
    // solve least squard problem for Ax=b as QR decomposition way if b is
    // [[b1],[b2],[b3]] form will return [[x1],[x2],[x3]] array form solution
    // else b is [b1,b2,b3] form will return [x1,x2,x3] array form solution
    function R_I(A) {
      A = jStat.copy(A);
      var size = A.length;
      var I = jStat.identity(size);
      jStat.arange(size - 1, -1, -1).forEach(function(i) {
        jStat.sliceAssign(
            I, { row: i }, jStat.divide(jStat.slice(I, { row: i }), A[i][i]));
        jStat.sliceAssign(
            A, { row: i }, jStat.divide(jStat.slice(A, { row: i }), A[i][i]));
        jStat.arange(i).forEach(function(j) {
          var c = jStat.multiply(A[j][i], -1);
          var Aj = jStat.slice(A, { row: j });
          var cAi = jStat.multiply(jStat.slice(A, { row: i }), c);
          jStat.sliceAssign(A, { row: j }, jStat.add(Aj, cAi));
          var Ij = jStat.slice(I, { row: j });
          var cIi = jStat.multiply(jStat.slice(I, { row: i }), c);
          jStat.sliceAssign(I, { row: j }, jStat.add(Ij, cIi));
        })
      });
      return I;
    }

    function qr_solve(A, b){
      var array_mode = false;
      if (b[0].length === undefined) {
        // [c1,c2,c3] mode
        b = b.map(function(x){ return [x] });
        array_mode = true;
      }
      var QR = jStat.QR(A);
      var Q = QR[0];
      var R = QR[1];
      var attrs = A[0].length;
      var Q1 = jStat.slice(Q,{col:{end:attrs}});
      var R1 = jStat.slice(R,{row:{end:attrs}});
      var RI = R_I(R1);
      var Q2 = jStat.transpose(Q1);

      if(Q2[0].length === undefined){
        Q2 = [Q2]; // The confusing jStat.multifly implementation threat nature process again.
      }

      var x = jStat.multiply(jStat.multiply(RI, Q2), b);

      if(x.length === undefined){
        x = [[x]]; // The confusing jStat.multifly implementation threat nature process again.
      }


      if (array_mode)
        return x.map(function(i){ return i[0] });
      return x;
    }

    return qr_solve;
  }()),

  jacobi: function jacobi(a) {
    var condition = 1;
    var n = a.length;
    var e = jStat.identity(n, n);
    var ev = [];
    var b, i, j, p, q, maxim, theta, s;
    // condition === 1 only if tolerance is not reached
    while (condition === 1) {
      maxim = a[0][1];
      p = 0;
      q = 1;
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
          if (i != j) {
            if (maxim < Math.abs(a[i][j])) {
              maxim = Math.abs(a[i][j]);
              p = i;
              q = j;
            }
          }
        }
      }
      if (a[p][p] === a[q][q])
        theta = (a[p][q] > 0) ? Math.PI / 4 : -Math.PI / 4;
      else
        theta = Math.atan(2 * a[p][q] / (a[p][p] - a[q][q])) / 2;
      s = jStat.identity(n, n);
      s[p][p] = Math.cos(theta);
      s[p][q] = -Math.sin(theta);
      s[q][p] = Math.sin(theta);
      s[q][q] = Math.cos(theta);
      // eigen vector matrix
      e = jStat.multiply(e, s);
      b = jStat.multiply(jStat.multiply(jStat.inv(s), a), s);
      a = b;
      condition = 0;
      for (i = 1; i < n; i++) {
        for (j = 1; j < n; j++) {
          if (i != j && Math.abs(a[i][j]) > 0.001) {
            condition = 1;
          }
        }
      }
    }
    for (i = 0; i < n; i++) ev.push(a[i][i]);
    //returns both the eigenvalue and eigenmatrix
    return [e, ev];
  },

  rungekutta: function rungekutta(f, h, p, t_j, u_j, order) {
    var k1, k2, u_j1, k3, k4;
    if (order === 2) {
      while (t_j <= p) {
        k1 = h * f(t_j, u_j);
        k2 = h * f(t_j + h, u_j + k1);
        u_j1 = u_j + (k1 + k2) / 2;
        u_j = u_j1;
        t_j = t_j + h;
      }
    }
    if (order === 4) {
      while (t_j <= p) {
        k1 = h * f(t_j, u_j);
        k2 = h * f(t_j + h / 2, u_j + k1 / 2);
        k3 = h * f(t_j + h / 2, u_j + k2 / 2);
        k4 = h * f(t_j +h, u_j + k3);
        u_j1 = u_j + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        u_j = u_j1;
        t_j = t_j + h;
      }
    }
    return u_j;
  },

  romberg: function romberg(f, a, b, order) {
    var i = 0;
    var h = (b - a) / 2;
    var x = [];
    var h1 = [];
    var g = [];
    var m, a1, j, k, I;
    while (i < order / 2) {
      I = f(a);
      for (j = a, k = 0; j <= b; j = j + h, k++) x[k] = j;
      m = x.length;
      for (j = 1; j < m - 1; j++) {
        I += (((j % 2) !== 0) ? 4 : 2) * f(x[j]);
      }
      I = (h / 3) * (I + f(b));
      g[i] = I;
      h /= 2;
      i++;
    }
    a1 = g.length;
    m = 1;
    while (a1 !== 1) {
      for (j = 0; j < a1 - 1; j++)
      h1[j] = ((Math.pow(4, m)) * g[j + 1] - g[j]) / (Math.pow(4, m) - 1);
      a1 = h1.length;
      g = h1;
      h1 = [];
      m++;
    }
    return g;
  },

  richardson: function richardson(X, f, x, h) {
    function pos(X, x) {
      var i = 0;
      var n = X.length;
      var p;
      for (; i < n; i++)
        if (X[i] === x) p = i;
      return p;
    }
    var h_min = Math.abs(x - X[pos(X, x) + 1]);
    var i = 0;
    var g = [];
    var h1 = [];
    var y1, y2, m, a, j;
    while (h >= h_min) {
      y1 = pos(X, x + h);
      y2 = pos(X, x);
      g[i] = (f[y1] - 2 * f[y2] + f[2 * y2 - y1]) / (h * h);
      h /= 2;
      i++;
    }
    a = g.length;
    m = 1;
    while (a != 1) {
      for (j = 0; j < a - 1; j++)
        h1[j] = ((Math.pow(4, m)) * g[j + 1] - g[j]) / (Math.pow(4, m) - 1);
      a = h1.length;
      g = h1;
      h1 = [];
      m++;
    }
    return g;
  },

  simpson: function simpson(f, a, b, n) {
    var h = (b - a) / n;
    var I = f(a);
    var x = [];
    var j = a;
    var k = 0;
    var i = 1;
    var m;
    for (; j <= b; j = j + h, k++)
      x[k] = j;
    m = x.length;
    for (; i < m - 1; i++) {
      I += ((i % 2 !== 0) ? 4 : 2) * f(x[i]);
    }
    return (h / 3) * (I + f(b));
  },

  hermite: function hermite(X, F, dF, value) {
    var n = X.length;
    var p = 0;
    var i = 0;
    var l = [];
    var dl = [];
    var A = [];
    var B = [];
    var j;
    for (; i < n; i++) {
      l[i] = 1;
      for (j = 0; j < n; j++) {
        if (i != j) l[i] *= (value - X[j]) / (X[i] - X[j]);
      }
      dl[i] = 0;
      for (j = 0; j < n; j++) {
        if (i != j) dl[i] += 1 / (X [i] - X[j]);
      }
      A[i] = (1 - 2 * (value - X[i]) * dl[i]) * (l[i] * l[i]);
      B[i] = (value - X[i]) * (l[i] * l[i]);
      p += (A[i] * F[i] + B[i] * dF[i]);
    }
    return p;
  },

  lagrange: function lagrange(X, F, value) {
    var p = 0;
    var i = 0;
    var j, l;
    var n = X.length;
    for (; i < n; i++) {
      l = F[i];
      for (j = 0; j < n; j++) {
        // calculating the lagrange polynomial L_i
        if (i != j) l *= (value - X[j]) / (X[i] - X[j]);
      }
      // adding the lagrange polynomials found above
      p += l;
    }
    return p;
  },

  cubic_spline: function cubic_spline(X, F, value) {
    var n = X.length;
    var i = 0, j;
    var A = [];
    var B = [];
    var alpha = [];
    var c = [];
    var h = [];
    var b = [];
    var d = [];
    for (; i < n - 1; i++)
      h[i] = X[i + 1] - X[i];
    alpha[0] = 0;
    for (i = 1; i < n - 1; i++) {
      alpha[i] = (3 / h[i]) * (F[i + 1] - F[i]) -
          (3 / h[i-1]) * (F[i] - F[i-1]);
    }
    for (i = 1; i < n - 1; i++) {
      A[i] = [];
      B[i] = [];
      A[i][i-1] = h[i-1];
      A[i][i] = 2 * (h[i - 1] + h[i]);
      A[i][i+1] = h[i];
      B[i][0] = alpha[i];
    }
    c = jStat.multiply(jStat.inv(A), B);
    for (j = 0; j < n - 1; j++) {
      b[j] = (F[j + 1] - F[j]) / h[j] - h[j] * (c[j + 1][0] + 2 * c[j][0]) / 3;
      d[j] = (c[j + 1][0] - c[j][0]) / (3 * h[j]);
    }
    for (j = 0; j < n; j++) {
      if (X[j] > value) break;
    }
    j -= 1;
    return F[j] + (value - X[j]) * b[j] + jStat.sq(value-X[j]) *
        c[j] + (value - X[j]) * jStat.sq(value - X[j]) * d[j];
  },

  gauss_quadrature: function gauss_quadrature() {
    throw new Error('gauss_quadrature not yet implemented');
  },

  PCA: function PCA(X) {
    var m = X.length;
    var n = X[0].length;
    var i = 0;
    var j, temp1;
    var u = [];
    var D = [];
    var result = [];
    var temp2 = [];
    var Y = [];
    var Bt = [];
    var B = [];
    var C = [];
    var V = [];
    var Vt = [];
    for (i = 0; i < m; i++) {
      u[i] = jStat.sum(X[i]) / n;
    }
    for (i = 0; i < n; i++) {
      B[i] = [];
      for(j = 0; j < m; j++) {
        B[i][j] = X[j][i] - u[j];
      }
    }
    B = jStat.transpose(B);
    for (i = 0; i < m; i++) {
      C[i] = [];
      for (j = 0; j < m; j++) {
        C[i][j] = (jStat.dot([B[i]], [B[j]])) / (n - 1);
      }
    }
    result = jStat.jacobi(C);
    V = result[0];
    D = result[1];
    Vt = jStat.transpose(V);
    for (i = 0; i < D.length; i++) {
      for (j = i; j < D.length; j++) {
        if(D[i] < D[j])  {
          temp1 = D[i];
          D[i] = D[j];
          D[j] = temp1;
          temp2 = Vt[i];
          Vt[i] = Vt[j];
          Vt[j] = temp2;
        }
      }
    }
    Bt = jStat.transpose(B);
    for (i = 0; i < m; i++) {
      Y[i] = [];
      for (j = 0; j < Bt.length; j++) {
        Y[i][j] = jStat.dot([Vt[i]], [Bt[j]]);
      }
    }
    return [X, D, Vt, Y];
  }
});

// extend jStat.fn with methods that require one argument
(function(funcs) {
  for (var i = 0; i < funcs.length; i++) (function(passfunc) {
    jStat.fn[passfunc] = function(arg, func) {
      var tmpthis = this;
      // check for callback
      if (func) {
        setTimeout(function() {
          func.call(tmpthis, jStat.fn[passfunc].call(tmpthis, arg));
        }, 15);
        return this;
      }
      if (typeof jStat[passfunc](this, arg) === 'number')
        return jStat[passfunc](this, arg);
      else
        return jStat(jStat[passfunc](this, arg));
    };
  }(funcs[i]));
}('add divide multiply subtract dot pow exp log abs norm angle'.split(' ')));

}(jStat, Math));
(function(jStat, Math) {

var slice = [].slice;
var isNumber = jStat.utils.isNumber;
var isArray = jStat.utils.isArray;

// flag==true denotes use of sample standard deviation
// Z Statistics
jStat.extend({
  // 2 different parameter lists:
  // (value, mean, sd)
  // (value, array, flag)
  zscore: function zscore() {
    var args = slice.call(arguments);
    if (isNumber(args[1])) {
      return (args[0] - args[1]) / args[2];
    }
    return (args[0] - jStat.mean(args[1])) / jStat.stdev(args[1], args[2]);
  },

  // 3 different paramter lists:
  // (value, mean, sd, sides)
  // (zscore, sides)
  // (value, array, sides, flag)
  ztest: function ztest() {
    var args = slice.call(arguments);
    var z;
    if (isArray(args[1])) {
      // (value, array, sides, flag)
      z = jStat.zscore(args[0],args[1],args[3]);
      return (args[2] === 1) ?
        (jStat.normal.cdf(-Math.abs(z), 0, 1)) :
        (jStat.normal.cdf(-Math.abs(z), 0, 1)*2);
    } else {
      if (args.length > 2) {
        // (value, mean, sd, sides)
        z = jStat.zscore(args[0],args[1],args[2]);
        return (args[3] === 1) ?
          (jStat.normal.cdf(-Math.abs(z),0,1)) :
          (jStat.normal.cdf(-Math.abs(z),0,1)* 2);
      } else {
        // (zscore, sides)
        z = args[0];
        return (args[1] === 1) ?
          (jStat.normal.cdf(-Math.abs(z),0,1)) :
          (jStat.normal.cdf(-Math.abs(z),0,1)*2);
      }
    }
  }
});

jStat.extend(jStat.fn, {
  zscore: function zscore(value, flag) {
    return (value - this.mean()) / this.stdev(flag);
  },

  ztest: function ztest(value, sides, flag) {
    var zscore = Math.abs(this.zscore(value, flag));
    return (sides === 1) ?
      (jStat.normal.cdf(-zscore, 0, 1)) :
      (jStat.normal.cdf(-zscore, 0, 1) * 2);
  }
});

// T Statistics
jStat.extend({
  // 2 parameter lists
  // (value, mean, sd, n)
  // (value, array)
  tscore: function tscore() {
    var args = slice.call(arguments);
    return (args.length === 4) ?
      ((args[0] - args[1]) / (args[2] / Math.sqrt(args[3]))) :
      ((args[0] - jStat.mean(args[1])) /
       (jStat.stdev(args[1], true) / Math.sqrt(args[1].length)));
  },

  // 3 different paramter lists:
  // (value, mean, sd, n, sides)
  // (tscore, n, sides)
  // (value, array, sides)
  ttest: function ttest() {
    var args = slice.call(arguments);
    var tscore;
    if (args.length === 5) {
      tscore = Math.abs(jStat.tscore(args[0], args[1], args[2], args[3]));
      return (args[4] === 1) ?
        (jStat.studentt.cdf(-tscore, args[3]-1)) :
        (jStat.studentt.cdf(-tscore, args[3]-1)*2);
    }
    if (isNumber(args[1])) {
      tscore = Math.abs(args[0])
      return (args[2] == 1) ?
        (jStat.studentt.cdf(-tscore, args[1]-1)) :
        (jStat.studentt.cdf(-tscore, args[1]-1) * 2);
    }
    tscore = Math.abs(jStat.tscore(args[0], args[1]))
    return (args[2] == 1) ?
      (jStat.studentt.cdf(-tscore, args[1].length-1)) :
      (jStat.studentt.cdf(-tscore, args[1].length-1) * 2);
  }
});

jStat.extend(jStat.fn, {
  tscore: function tscore(value) {
    return (value - this.mean()) / (this.stdev(true) / Math.sqrt(this.cols()));
  },

  ttest: function ttest(value, sides) {
    return (sides === 1) ?
      (1 - jStat.studentt.cdf(Math.abs(this.tscore(value)), this.cols()-1)) :
      (jStat.studentt.cdf(-Math.abs(this.tscore(value)), this.cols()-1)*2);
  }
});

// F Statistics
jStat.extend({
  // Paramter list is as follows:
  // (array1, array2, array3, ...)
  // or it is an array of arrays
  // array of arrays conversion
  anovafscore: function anovafscore() {
    var args = slice.call(arguments),
    expVar, sample, sampMean, sampSampMean, tmpargs, unexpVar, i, j;
    if (args.length === 1) {
      tmpargs = new Array(args[0].length);
      for (i = 0; i < args[0].length; i++) {
        tmpargs[i] = args[0][i];
      }
      args = tmpargs;
    }
    // Builds sample array
    sample = new Array();
    for (i = 0; i < args.length; i++) {
      sample = sample.concat(args[i]);
    }
    sampMean = jStat.mean(sample);
    // Computes the explained variance
    expVar = 0;
    for (i = 0; i < args.length; i++) {
      expVar = expVar + args[i].length * Math.pow(jStat.mean(args[i]) - sampMean, 2);
    }
    expVar /= (args.length - 1);
    // Computes unexplained variance
    unexpVar = 0;
    for (i = 0; i < args.length; i++) {
      sampSampMean = jStat.mean(args[i]);
      for (j = 0; j < args[i].length; j++) {
        unexpVar += Math.pow(args[i][j] - sampSampMean, 2);
      }
    }
    unexpVar /= (sample.length - args.length);
    return expVar / unexpVar;
  },

  // 2 different paramter setups
  // (array1, array2, array3, ...)
  // (anovafscore, df1, df2)
  anovaftest: function anovaftest() {
    var args = slice.call(arguments),
    df1, df2, n, i;
    if (isNumber(args[0])) {
      return 1 - jStat.centralF.cdf(args[0], args[1], args[2]);
    }
    var anovafscore = jStat.anovafscore(args);
    df1 = args.length - 1;
    n = 0;
    for (i = 0; i < args.length; i++) {
      n = n + args[i].length;
    }
    df2 = n - df1 - 1;
    return 1 - jStat.centralF.cdf(anovafscore, df1, df2);
  },

  ftest: function ftest(fscore, df1, df2) {
    return 1 - jStat.centralF.cdf(fscore, df1, df2);
  }
});

jStat.extend(jStat.fn, {
  anovafscore: function anovafscore() {
    return jStat.anovafscore(this.toArray());
  },

  anovaftes: function anovaftes() {
    var n = 0;
    var i;
    for (i = 0; i < this.length; i++) {
      n = n + this[i].length;
    }
    return jStat.ftest(this.anovafscore(), this.length - 1, n - this.length);
  }
});

// Tukey's range test
jStat.extend({
  // 2 parameter lists
  // (mean1, mean2, n1, n2, sd)
  // (array1, array2, sd)
  qscore: function qscore() {
    var args = slice.call(arguments);
    var mean1, mean2, n1, n2, sd;
    if (isNumber(args[0])) {
        mean1 = args[0];
        mean2 = args[1];
        n1 = args[2];
        n2 = args[3];
        sd = args[4];
    } else {
        mean1 = jStat.mean(args[0]);
        mean2 = jStat.mean(args[1]);
        n1 = args[0].length;
        n2 = args[1].length;
        sd = args[2];
    }
    return Math.abs(mean1 - mean2) / (sd * Math.sqrt((1 / n1 + 1 / n2) / 2));
  },

  // 3 different parameter lists:
  // (qscore, n, k)
  // (mean1, mean2, n1, n2, sd, n, k)
  // (array1, array2, sd, n, k)
  qtest: function qtest() {
    var args = slice.call(arguments);

    var qscore;
    if (args.length === 3) {
      qscore = args[0];
      args = args.slice(1);
    } else if (args.length === 7) {
      qscore = jStat.qscore(args[0], args[1], args[2], args[3], args[4]);
      args = args.slice(5);
    } else {
      qscore = jStat.qscore(args[0], args[1], args[2]);
      args = args.slice(3);
    }

    var n = args[0];
    var k = args[1];

    return 1 - jStat.tukey.cdf(qscore, k, n - k);
  },

  tukeyhsd: function tukeyhsd(arrays) {
    var sd = jStat.pooledstdev(arrays);
    var means = arrays.map(function (arr) {return jStat.mean(arr);});
    var n = arrays.reduce(function (n, arr) {return n + arr.length;}, 0);

    var results = [];
    for (var i = 0; i < arrays.length; ++i) {
        for (var j = i + 1; j < arrays.length; ++j) {
            var p = jStat.qtest(means[i], means[j], arrays[i].length, arrays[j].length, sd, n, arrays.length);
            results.push([[i, j], p]);
        }
    }

    return results;
  }
});

// Error Bounds
jStat.extend({
  // 2 different parameter setups
  // (value, alpha, sd, n)
  // (value, alpha, array)
  normalci: function normalci() {
    var args = slice.call(arguments),
    ans = new Array(2),
    change;
    if (args.length === 4) {
      change = Math.abs(jStat.normal.inv(args[1] / 2, 0, 1) *
                        args[2] / Math.sqrt(args[3]));
    } else {
      change = Math.abs(jStat.normal.inv(args[1] / 2, 0, 1) *
                        jStat.stdev(args[2]) / Math.sqrt(args[2].length));
    }
    ans[0] = args[0] - change;
    ans[1] = args[0] + change;
    return ans;
  },

  // 2 different parameter setups
  // (value, alpha, sd, n)
  // (value, alpha, array)
  tci: function tci() {
    var args = slice.call(arguments),
    ans = new Array(2),
    change;
    if (args.length === 4) {
      change = Math.abs(jStat.studentt.inv(args[1] / 2, args[3] - 1) *
                        args[2] / Math.sqrt(args[3]));
    } else {
      change = Math.abs(jStat.studentt.inv(args[1] / 2, args[2].length - 1) *
                        jStat.stdev(args[2], true) / Math.sqrt(args[2].length));
    }
    ans[0] = args[0] - change;
    ans[1] = args[0] + change;
    return ans;
  },

  significant: function significant(pvalue, alpha) {
    return pvalue < alpha;
  }
});

jStat.extend(jStat.fn, {
  normalci: function normalci(value, alpha) {
    return jStat.normalci(value, alpha, this.toArray());
  },

  tci: function tci(value, alpha) {
    return jStat.tci(value, alpha, this.toArray());
  }
});

// internal method for calculating the z-score for a difference of proportions test
function differenceOfProportions(p1, n1, p2, n2) {
  if (p1 > 1 || p2 > 1 || p1 <= 0 || p2 <= 0) {
    throw new Error("Proportions should be greater than 0 and less than 1")
  }
  var pooled = (p1 * n1 + p2 * n2) / (n1 + n2);
  var se = Math.sqrt(pooled * (1 - pooled) * ((1/n1) + (1/n2)));
  return (p1 - p2) / se;
}

// Difference of Proportions
jStat.extend(jStat.fn, {
  oneSidedDifferenceOfProportions: function oneSidedDifferenceOfProportions(p1, n1, p2, n2) {
    var z = differenceOfProportions(p1, n1, p2, n2);
    return jStat.ztest(z, 1);
  },

  twoSidedDifferenceOfProportions: function twoSidedDifferenceOfProportions(p1, n1, p2, n2) {
    var z = differenceOfProportions(p1, n1, p2, n2);
    return jStat.ztest(z, 2);
  }
});

}(jStat, Math));
jStat.models = (function(){
  function sub_regress(exog) {
    var var_count = exog[0].length;
    var modelList = jStat.arange(var_count).map(function(endog_index) {
      var exog_index =
          jStat.arange(var_count).filter(function(i){return i!==endog_index});
      return ols(jStat.col(exog, endog_index).map(function(x){ return x[0] }),
                 jStat.col(exog, exog_index))
    });
    return modelList;
  }

  // do OLS model regress
  // exog have include const columns ,it will not generate it .In fact, exog is
  // "design matrix" look at
  //https://en.wikipedia.org/wiki/Design_matrix
  function ols(endog, exog) {
    var nobs = endog.length;
    var df_model = exog[0].length - 1;
    var df_resid = nobs-df_model - 1;
    var coef = jStat.lstsq(exog, endog);
    var predict =
        jStat.multiply(exog, coef.map(function(x) { return [x] }))
            .map(function(p) { return p[0] });
    var resid = jStat.subtract(endog, predict);
    var ybar = jStat.mean(endog);
    // constant cause problem
    // var SST = jStat.sum(endog.map(function(y) {
    //   return Math.pow(y-ybar,2);
    // }));
    var SSE = jStat.sum(predict.map(function(f) {
      return Math.pow(f - ybar, 2);
    }));
    var SSR = jStat.sum(endog.map(function(y, i) {
      return Math.pow(y - predict[i], 2);
    }));
    var SST = SSE + SSR;
    var R2 = (SSE / SST);
    return {
        exog:exog,
        endog:endog,
        nobs:nobs,
        df_model:df_model,
        df_resid:df_resid,
        coef:coef,
        predict:predict,
        resid:resid,
        ybar:ybar,
        SST:SST,
        SSE:SSE,
        SSR:SSR,
        R2:R2
    };
  }

  // H0: b_I=0
  // H1: b_I!=0
  function t_test(model) {
    var subModelList = sub_regress(model.exog);
    //var sigmaHat=jStat.stdev(model.resid);
    var sigmaHat = Math.sqrt(model.SSR / (model.df_resid));
    var seBetaHat = subModelList.map(function(mod) {
      var SST = mod.SST;
      var R2 = mod.R2;
      return sigmaHat / Math.sqrt(SST * (1 - R2));
    });
    var tStatistic = model.coef.map(function(coef, i) {
      return (coef - 0) / seBetaHat[i];
    });
    var pValue = tStatistic.map(function(t) {
      var leftppf = jStat.studentt.cdf(t, model.df_resid);
      return (leftppf > 0.5 ? 1 - leftppf : leftppf) * 2;
    });
    var c = jStat.studentt.inv(0.975, model.df_resid);
    var interval95 = model.coef.map(function(coef, i) {
      var d = c * seBetaHat[i];
      return [coef - d, coef + d];
    })
    return {
        se: seBetaHat,
        t: tStatistic,
        p: pValue,
        sigmaHat: sigmaHat,
        interval95: interval95
    };
  }

  function F_test(model) {
    var F_statistic =
        (model.R2 / model.df_model) / ((1 - model.R2) / model.df_resid);
    var fcdf = function(x, n1, n2) {
      return jStat.beta.cdf(x / (n2 / n1 + x), n1 / 2, n2 / 2)
    }
    var pvalue = 1 - fcdf(F_statistic, model.df_model, model.df_resid);
    return { F_statistic: F_statistic, pvalue: pvalue };
  }

  function ols_wrap(endog, exog) {
    var model = ols(endog,exog);
    var ttest = t_test(model);
    var ftest = F_test(model);
    // Provide the Wherry / Ezekiel / McNemar / Cohen Adjusted R^2
    // Which matches the 'adjusted R^2' provided by R's lm package
    var adjust_R2 =
        1 - (1 - model.R2) * ((model.nobs - 1) / (model.df_resid));
    model.t = ttest;
    model.f = ftest;
    model.adjust_R2 = adjust_R2;
    return model;
  }

  return { ols: ols_wrap };
})();
//To regress, simply build X matrix
//(append column of 1's) using
//buildxmatrix and build the Y
//matrix using buildymatrix
//(simply the transpose)
//and run regress.



//Regressions

jStat.extend({
  buildxmatrix: function buildxmatrix(){
    //Parameters will be passed in as such
    //(array1,array2,array3,...)
    //as (x1,x2,x3,...)
    //needs to be (1,x1,x2,x3,...)
    var matrixRows = new Array(arguments.length);
    for(var i=0;i<arguments.length;i++){
      var array = [1];
      matrixRows[i]= array.concat(arguments[i]);
    }
    return jStat(matrixRows);

  },

  builddxmatrix: function builddxmatrix() {
    //Paramters will be passed in as such
    //([array1,array2,...]
    var matrixRows = new Array(arguments[0].length);
    for(var i=0;i<arguments[0].length;i++){
      var array = [1]
      matrixRows[i]= array.concat(arguments[0][i]);
    }
    return jStat(matrixRows);

  },

  buildjxmatrix: function buildjxmatrix(jMat) {
    //Builds from jStat Matrix
    var pass = new Array(jMat.length)
    for(var i=0;i<jMat.length;i++){
      pass[i] = jMat[i];
    }
    return jStat.builddxmatrix(pass);

  },

  buildymatrix: function buildymatrix(array){
    return jStat(array).transpose();
  },

  buildjymatrix: function buildjymatrix(jMat){
    return jMat.transpose();
  },

  matrixmult: function matrixmult(A,B){
    var i, j, k, result, sum;
    if (A.cols() == B.rows()) {
      if(B.rows()>1){
        result = [];
        for (i = 0; i < A.rows(); i++) {
          result[i] = [];
          for (j = 0; j < B.cols(); j++) {
            sum = 0;
            for (k = 0; k < A.cols(); k++) {
              sum += A.toArray()[i][k] * B.toArray()[k][j];
            }
            result[i][j] = sum;
          }
        }
        return jStat(result);
      }
      result = [];
      for (i = 0; i < A.rows(); i++) {
        result[i] = [];
        for (j = 0; j < B.cols(); j++) {
          sum = 0;
          for (k = 0; k < A.cols(); k++) {
            sum += A.toArray()[i][k] * B.toArray()[j];
          }
          result[i][j] = sum;
        }
      }
      return jStat(result);
    }
  },

  //regress and regresst to be fixed

  regress: function regress(jMatX,jMatY){
    //print("regressin!");
    //print(jMatX.toArray());
    var innerinv = jStat.xtranspxinv(jMatX);
    //print(innerinv);
    var xtransp = jMatX.transpose();
    var next = jStat.matrixmult(jStat(innerinv),xtransp);
    return jStat.matrixmult(next,jMatY);

  },

  regresst: function regresst(jMatX,jMatY,sides){
    var beta = jStat.regress(jMatX,jMatY);

    var compile = {};
    compile.anova = {};
    var jMatYBar = jStat.jMatYBar(jMatX, beta);
    compile.yBar = jMatYBar;
    var yAverage = jMatY.mean();
    compile.anova.residuals = jStat.residuals(jMatY, jMatYBar);

    compile.anova.ssr = jStat.ssr(jMatYBar, yAverage);
    compile.anova.msr = compile.anova.ssr / (jMatX[0].length - 1);

    compile.anova.sse = jStat.sse(jMatY, jMatYBar);
    compile.anova.mse =
        compile.anova.sse / (jMatY.length - (jMatX[0].length - 1) - 1);

    compile.anova.sst = jStat.sst(jMatY, yAverage);
    compile.anova.mst = compile.anova.sst / (jMatY.length - 1);

    compile.anova.r2 = 1 - (compile.anova.sse / compile.anova.sst);
    if (compile.anova.r2 < 0) compile.anova.r2 = 0;

    compile.anova.fratio = compile.anova.msr / compile.anova.mse;
    compile.anova.pvalue =
        jStat.anovaftest(compile.anova.fratio,
                         jMatX[0].length - 1,
                         jMatY.length - (jMatX[0].length - 1) - 1);

    compile.anova.rmse = Math.sqrt(compile.anova.mse);

    compile.anova.r2adj = 1 - (compile.anova.mse / compile.anova.mst);
    if (compile.anova.r2adj < 0) compile.anova.r2adj = 0;

    compile.stats = new Array(jMatX[0].length);
    var covar = jStat.xtranspxinv(jMatX);
    var sds, ts, ps;

    for(var i=0; i<beta.length;i++){
      sds=Math.sqrt(compile.anova.mse * Math.abs(covar[i][i]));
      ts= Math.abs(beta[i] / sds);
      ps= jStat.ttest(ts, jMatY.length - jMatX[0].length - 1, sides);

      compile.stats[i]=[beta[i], sds, ts, ps];
    }

    compile.regress = beta;
    return compile;
  },

  xtranspx: function xtranspx(jMatX){
    return jStat.matrixmult(jMatX.transpose(),jMatX);
  },


  xtranspxinv: function xtranspxinv(jMatX){
    var inner = jStat.matrixmult(jMatX.transpose(),jMatX);
    var innerinv = jStat.inv(inner);
    return innerinv;
  },

  jMatYBar: function jMatYBar(jMatX, beta) {
    var yBar = jStat.matrixmult(jMatX, beta);
    return new jStat(yBar);
  },

  residuals: function residuals(jMatY, jMatYBar) {
    return jStat.matrixsubtract(jMatY, jMatYBar);
  },

  ssr: function ssr(jMatYBar, yAverage) {
    var ssr = 0;
    for(var i = 0; i < jMatYBar.length; i++) {
      ssr += Math.pow(jMatYBar[i] - yAverage, 2);
    }
    return ssr;
  },

  sse: function sse(jMatY, jMatYBar) {
    var sse = 0;
    for(var i = 0; i < jMatY.length; i++) {
      sse += Math.pow(jMatY[i] - jMatYBar[i], 2);
    }
    return sse;
  },

  sst: function sst(jMatY, yAverage) {
    var sst = 0;
    for(var i = 0; i < jMatY.length; i++) {
      sst += Math.pow(jMatY[i] - yAverage, 2);
    }
    return sst;
  },

  matrixsubtract: function matrixsubtract(A,B){
    var ans = new Array(A.length);
    for(var i=0;i<A.length;i++){
      ans[i] = new Array(A[i].length);
      for(var j=0;j<A[i].length;j++){
        ans[i][j]=A[i][j]-B[i][j];
      }
    }
    return jStat(ans);
  }
});
  // Make it compatible with previous version.
  jStat.jStat = jStat;

  return jStat;
});


/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	/* webpack/runtime/define property getters */
/******/ 	!function() {
/******/ 		// define getter functions for harmony exports
/******/ 		__webpack_require__.d = function(exports, definition) {
/******/ 			for(var key in definition) {
/******/ 				if(__webpack_require__.o(definition, key) && !__webpack_require__.o(exports, key)) {
/******/ 					Object.defineProperty(exports, key, { enumerable: true, get: definition[key] });
/******/ 				}
/******/ 			}
/******/ 		};
/******/ 	}();
/******/ 	
/******/ 	/* webpack/runtime/hasOwnProperty shorthand */
/******/ 	!function() {
/******/ 		__webpack_require__.o = function(obj, prop) { return Object.prototype.hasOwnProperty.call(obj, prop); }
/******/ 	}();
/******/ 	
/******/ 	/* webpack/runtime/make namespace object */
/******/ 	!function() {
/******/ 		// define __esModule on exports
/******/ 		__webpack_require__.r = function(exports) {
/******/ 			if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 				Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 			}
/******/ 			Object.defineProperty(exports, '__esModule', { value: true });
/******/ 		};
/******/ 	}();
/******/ 	
/************************************************************************/
var __webpack_exports__ = {};
// This entry need to be wrapped in an IIFE because it need to be in strict mode.
!function() {
"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  "default": function() { return /* binding */ src_formula; }
});

// NAMESPACE OBJECT: ./node_modules/@formulajs/formulajs/lib/esm/index.mjs
var esm_namespaceObject = {};
__webpack_require__.r(esm_namespaceObject);
__webpack_require__.d(esm_namespaceObject, {
  "ABS": function() { return ABS; },
  "ACCRINT": function() { return ACCRINT; },
  "ACCRINTM": function() { return ACCRINTM; },
  "ACOS": function() { return ACOS; },
  "ACOSH": function() { return ACOSH; },
  "ACOT": function() { return ACOT; },
  "ACOTH": function() { return ACOTH; },
  "AGGREGATE": function() { return AGGREGATE; },
  "AMORDEGRC": function() { return AMORDEGRC; },
  "AMORLINC": function() { return AMORLINC; },
  "AND": function() { return AND; },
  "ARABIC": function() { return ARABIC; },
  "ASC": function() { return ASC; },
  "ASIN": function() { return ASIN; },
  "ASINH": function() { return ASINH; },
  "ATAN": function() { return ATAN; },
  "ATAN2": function() { return ATAN2; },
  "ATANH": function() { return ATANH; },
  "AVEDEV": function() { return AVEDEV; },
  "AVERAGE": function() { return AVERAGE; },
  "AVERAGEA": function() { return AVERAGEA; },
  "AVERAGEIF": function() { return AVERAGEIF; },
  "AVERAGEIFS": function() { return AVERAGEIFS; },
  "BAHTTEXT": function() { return BAHTTEXT; },
  "BASE": function() { return BASE; },
  "BESSELI": function() { return BESSELI; },
  "BESSELJ": function() { return BESSELJ; },
  "BESSELK": function() { return BESSELK; },
  "BESSELY": function() { return BESSELY; },
  "BETA": function() { return BETA; },
  "BETADIST": function() { return BETADIST; },
  "BETAINV": function() { return BETAINV; },
  "BIN2DEC": function() { return BIN2DEC; },
  "BIN2HEX": function() { return BIN2HEX; },
  "BIN2OCT": function() { return BIN2OCT; },
  "BINOM": function() { return BINOM; },
  "BINOMDIST": function() { return BINOMDIST; },
  "BITAND": function() { return BITAND; },
  "BITLSHIFT": function() { return BITLSHIFT; },
  "BITOR": function() { return BITOR; },
  "BITRSHIFT": function() { return BITRSHIFT; },
  "BITXOR": function() { return BITXOR; },
  "CEILING": function() { return CEILING; },
  "CEILINGMATH": function() { return CEILINGMATH; },
  "CEILINGPRECISE": function() { return CEILINGPRECISE; },
  "CELL": function() { return CELL; },
  "CHAR": function() { return CHAR; },
  "CHIDIST": function() { return CHIDIST; },
  "CHIDISTRT": function() { return CHIDISTRT; },
  "CHIINV": function() { return CHIINV; },
  "CHIINVRT": function() { return CHIINVRT; },
  "CHISQ": function() { return CHISQ; },
  "CHITEST": function() { return CHITEST; },
  "CHOOSE": function() { return CHOOSE; },
  "CLEAN": function() { return CLEAN; },
  "CODE": function() { return CODE; },
  "COLUMN": function() { return COLUMN; },
  "COLUMNS": function() { return COLUMNS; },
  "COMBIN": function() { return COMBIN; },
  "COMBINA": function() { return COMBINA; },
  "COMPLEX": function() { return COMPLEX; },
  "CONCAT": function() { return CONCAT; },
  "CONCATENATE": function() { return CONCATENATE; },
  "CONFIDENCE": function() { return CONFIDENCE; },
  "CONVERT": function() { return CONVERT; },
  "CORREL": function() { return CORREL; },
  "COS": function() { return COS; },
  "COSH": function() { return COSH; },
  "COT": function() { return COT; },
  "COTH": function() { return COTH; },
  "COUNT": function() { return COUNT; },
  "COUNTA": function() { return COUNTA; },
  "COUNTBLANK": function() { return COUNTBLANK; },
  "COUNTIF": function() { return COUNTIF; },
  "COUNTIFS": function() { return COUNTIFS; },
  "COUPDAYBS": function() { return COUPDAYBS; },
  "COUPDAYS": function() { return COUPDAYS; },
  "COUPDAYSNC": function() { return COUPDAYSNC; },
  "COUPNCD": function() { return COUPNCD; },
  "COUPNUM": function() { return COUPNUM; },
  "COUPPCD": function() { return COUPPCD; },
  "COVAR": function() { return COVAR; },
  "COVARIANCE": function() { return COVARIANCE; },
  "COVARIANCEP": function() { return COVARIANCEP; },
  "COVARIANCES": function() { return COVARIANCES; },
  "CRITBINOM": function() { return CRITBINOM; },
  "CSC": function() { return CSC; },
  "CSCH": function() { return CSCH; },
  "CUMIPMT": function() { return CUMIPMT; },
  "CUMPRINC": function() { return CUMPRINC; },
  "DATE": function() { return DATE; },
  "DATEDIF": function() { return DATEDIF; },
  "DATEVALUE": function() { return DATEVALUE; },
  "DAVERAGE": function() { return DAVERAGE; },
  "DAY": function() { return DAY; },
  "DAYS": function() { return DAYS; },
  "DAYS360": function() { return DAYS360; },
  "DB": function() { return DB; },
  "DBCS": function() { return DBCS; },
  "DCOUNT": function() { return DCOUNT; },
  "DCOUNTA": function() { return DCOUNTA; },
  "DDB": function() { return DDB; },
  "DEC2BIN": function() { return DEC2BIN; },
  "DEC2HEX": function() { return DEC2HEX; },
  "DEC2OCT": function() { return DEC2OCT; },
  "DECIMAL": function() { return DECIMAL; },
  "DEGREES": function() { return DEGREES; },
  "DELTA": function() { return DELTA; },
  "DEVSQ": function() { return DEVSQ; },
  "DGET": function() { return DGET; },
  "DISC": function() { return DISC; },
  "DMAX": function() { return DMAX; },
  "DMIN": function() { return DMIN; },
  "DOLLAR": function() { return DOLLAR; },
  "DOLLARDE": function() { return DOLLARDE; },
  "DOLLARFR": function() { return DOLLARFR; },
  "DPRODUCT": function() { return DPRODUCT; },
  "DSTDEV": function() { return DSTDEV; },
  "DSTDEVP": function() { return DSTDEVP; },
  "DSUM": function() { return DSUM; },
  "DURATION": function() { return DURATION; },
  "DVAR": function() { return DVAR; },
  "DVARP": function() { return DVARP; },
  "EDATE": function() { return EDATE; },
  "EFFECT": function() { return EFFECT; },
  "EOMONTH": function() { return EOMONTH; },
  "ERF": function() { return ERF; },
  "ERFC": function() { return ERFC; },
  "ERFCPRECISE": function() { return ERFCPRECISE; },
  "ERFPRECISE": function() { return ERFPRECISE; },
  "ERROR": function() { return ERROR; },
  "EVEN": function() { return EVEN; },
  "EXACT": function() { return EXACT; },
  "EXP": function() { return EXP; },
  "EXPON": function() { return EXPON; },
  "EXPONDIST": function() { return EXPONDIST; },
  "F": function() { return F; },
  "FACT": function() { return FACT; },
  "FACTDOUBLE": function() { return FACTDOUBLE; },
  "FALSE": function() { return FALSE; },
  "FDIST": function() { return FDIST; },
  "FDISTRT": function() { return FDISTRT; },
  "FIND": function() { return FIND; },
  "FINV": function() { return FINV; },
  "FINVRT": function() { return FINVRT; },
  "FISHER": function() { return FISHER; },
  "FISHERINV": function() { return FISHERINV; },
  "FIXED": function() { return FIXED; },
  "FLOOR": function() { return FLOOR; },
  "FLOORMATH": function() { return FLOORMATH; },
  "FLOORPRECISE": function() { return FLOORPRECISE; },
  "FORECAST": function() { return FORECAST; },
  "FREQUENCY": function() { return FREQUENCY; },
  "FTEST": function() { return FTEST; },
  "FV": function() { return FV; },
  "FVSCHEDULE": function() { return FVSCHEDULE; },
  "GAMMA": function() { return GAMMA; },
  "GAMMADIST": function() { return GAMMADIST; },
  "GAMMAINV": function() { return GAMMAINV; },
  "GAMMALN": function() { return GAMMALN; },
  "GAMMALNPRECISE": function() { return GAMMALNPRECISE; },
  "GAUSS": function() { return GAUSS; },
  "GCD": function() { return GCD; },
  "GEOMEAN": function() { return GEOMEAN; },
  "GESTEP": function() { return GESTEP; },
  "GROWTH": function() { return GROWTH; },
  "HARMEAN": function() { return HARMEAN; },
  "HEX2BIN": function() { return HEX2BIN; },
  "HEX2DEC": function() { return HEX2DEC; },
  "HEX2OCT": function() { return HEX2OCT; },
  "HLOOKUP": function() { return HLOOKUP; },
  "HOUR": function() { return HOUR; },
  "HYPGEOM": function() { return HYPGEOM; },
  "HYPGEOMDIST": function() { return HYPGEOMDIST; },
  "IF": function() { return IF; },
  "IFERROR": function() { return IFERROR; },
  "IFNA": function() { return IFNA; },
  "IFS": function() { return IFS; },
  "IMABS": function() { return IMABS; },
  "IMAGINARY": function() { return IMAGINARY; },
  "IMARGUMENT": function() { return IMARGUMENT; },
  "IMCONJUGATE": function() { return IMCONJUGATE; },
  "IMCOS": function() { return IMCOS; },
  "IMCOSH": function() { return IMCOSH; },
  "IMCOT": function() { return IMCOT; },
  "IMCSC": function() { return IMCSC; },
  "IMCSCH": function() { return IMCSCH; },
  "IMDIV": function() { return IMDIV; },
  "IMEXP": function() { return IMEXP; },
  "IMLN": function() { return IMLN; },
  "IMLOG10": function() { return IMLOG10; },
  "IMLOG2": function() { return IMLOG2; },
  "IMPOWER": function() { return IMPOWER; },
  "IMPRODUCT": function() { return IMPRODUCT; },
  "IMREAL": function() { return IMREAL; },
  "IMSEC": function() { return IMSEC; },
  "IMSECH": function() { return IMSECH; },
  "IMSIN": function() { return IMSIN; },
  "IMSINH": function() { return IMSINH; },
  "IMSQRT": function() { return IMSQRT; },
  "IMSUB": function() { return IMSUB; },
  "IMSUM": function() { return IMSUM; },
  "IMTAN": function() { return IMTAN; },
  "INDEX": function() { return INDEX; },
  "INFO": function() { return INFO; },
  "INT": function() { return INT; },
  "INTERCEPT": function() { return INTERCEPT; },
  "INTRATE": function() { return INTRATE; },
  "IPMT": function() { return IPMT; },
  "IRR": function() { return IRR; },
  "ISBLANK": function() { return ISBLANK; },
  "ISERR": function() { return ISERR; },
  "ISERROR": function() { return ISERROR; },
  "ISEVEN": function() { return ISEVEN; },
  "ISFORMULA": function() { return ISFORMULA; },
  "ISLOGICAL": function() { return ISLOGICAL; },
  "ISNA": function() { return ISNA; },
  "ISNONTEXT": function() { return ISNONTEXT; },
  "ISNUMBER": function() { return ISNUMBER; },
  "ISO": function() { return ISO; },
  "ISODD": function() { return ISODD; },
  "ISOWEEKNUM": function() { return ISOWEEKNUM; },
  "ISPMT": function() { return ISPMT; },
  "ISREF": function() { return ISREF; },
  "ISTEXT": function() { return ISTEXT; },
  "KURT": function() { return KURT; },
  "LARGE": function() { return LARGE; },
  "LCM": function() { return LCM; },
  "LEFT": function() { return LEFT; },
  "LEN": function() { return LEN; },
  "LINEST": function() { return LINEST; },
  "LN": function() { return LN; },
  "LOG": function() { return LOG; },
  "LOG10": function() { return LOG10; },
  "LOGEST": function() { return LOGEST; },
  "LOGINV": function() { return LOGINV; },
  "LOGNORM": function() { return LOGNORM; },
  "LOGNORMDIST": function() { return LOGNORMDIST; },
  "LOGNORMINV": function() { return LOGNORMINV; },
  "LOOKUP": function() { return LOOKUP; },
  "LOWER": function() { return LOWER; },
  "MATCH": function() { return MATCH; },
  "MAX": function() { return MAX; },
  "MAXA": function() { return MAXA; },
  "MDURATION": function() { return MDURATION; },
  "MEDIAN": function() { return MEDIAN; },
  "MID": function() { return MID; },
  "MIN": function() { return MIN; },
  "MINA": function() { return MINA; },
  "MINUTE": function() { return MINUTE; },
  "MIRR": function() { return MIRR; },
  "MMULT": function() { return MMULT; },
  "MOD": function() { return MOD; },
  "MODE": function() { return MODE; },
  "MODEMULT": function() { return MODEMULT; },
  "MODESNGL": function() { return MODESNGL; },
  "MONTH": function() { return MONTH; },
  "MROUND": function() { return MROUND; },
  "MULTINOMIAL": function() { return MULTINOMIAL; },
  "MUNIT": function() { return MUNIT; },
  "N": function() { return N; },
  "NA": function() { return NA; },
  "NEGBINOM": function() { return NEGBINOM; },
  "NEGBINOMDIST": function() { return NEGBINOMDIST; },
  "NETWORKDAYS": function() { return NETWORKDAYS; },
  "NETWORKDAYSINTL": function() { return NETWORKDAYSINTL; },
  "NOMINAL": function() { return NOMINAL; },
  "NORM": function() { return NORM; },
  "NORMDIST": function() { return NORMDIST; },
  "NORMINV": function() { return NORMINV; },
  "NORMSDIST": function() { return NORMSDIST; },
  "NORMSINV": function() { return NORMSINV; },
  "NOT": function() { return NOT; },
  "NOW": function() { return NOW; },
  "NPER": function() { return NPER; },
  "NPV": function() { return NPV; },
  "NUMBERVALUE": function() { return NUMBERVALUE; },
  "OCT2BIN": function() { return OCT2BIN; },
  "OCT2DEC": function() { return OCT2DEC; },
  "OCT2HEX": function() { return OCT2HEX; },
  "ODD": function() { return ODD; },
  "ODDFPRICE": function() { return ODDFPRICE; },
  "ODDFYIELD": function() { return ODDFYIELD; },
  "ODDLPRICE": function() { return ODDLPRICE; },
  "ODDLYIELD": function() { return ODDLYIELD; },
  "OR": function() { return OR; },
  "PDURATION": function() { return PDURATION; },
  "PEARSON": function() { return PEARSON; },
  "PERCENTILE": function() { return PERCENTILE; },
  "PERCENTILEEXC": function() { return PERCENTILEEXC; },
  "PERCENTILEINC": function() { return PERCENTILEINC; },
  "PERCENTRANK": function() { return PERCENTRANK; },
  "PERCENTRANKEXC": function() { return PERCENTRANKEXC; },
  "PERCENTRANKINC": function() { return PERCENTRANKINC; },
  "PERMUT": function() { return PERMUT; },
  "PERMUTATIONA": function() { return PERMUTATIONA; },
  "PHI": function() { return PHI; },
  "PI": function() { return PI; },
  "PMT": function() { return PMT; },
  "POISSON": function() { return POISSON; },
  "POISSONDIST": function() { return POISSONDIST; },
  "POWER": function() { return POWER; },
  "PPMT": function() { return PPMT; },
  "PRICE": function() { return PRICE; },
  "PRICEDISC": function() { return PRICEDISC; },
  "PRICEMAT": function() { return PRICEMAT; },
  "PROB": function() { return PROB; },
  "PRODUCT": function() { return PRODUCT; },
  "PRONETIC": function() { return PRONETIC; },
  "PROPER": function() { return PROPER; },
  "PV": function() { return PV; },
  "QUARTILE": function() { return QUARTILE; },
  "QUARTILEEXC": function() { return QUARTILEEXC; },
  "QUARTILEINC": function() { return QUARTILEINC; },
  "QUOTIENT": function() { return QUOTIENT; },
  "RADIANS": function() { return RADIANS; },
  "RAND": function() { return RAND; },
  "RANDBETWEEN": function() { return RANDBETWEEN; },
  "RANK": function() { return RANK; },
  "RANKAVG": function() { return RANKAVG; },
  "RANKEQ": function() { return RANKEQ; },
  "RATE": function() { return RATE; },
  "RECEIVED": function() { return RECEIVED; },
  "REPLACE": function() { return REPLACE; },
  "REPT": function() { return REPT; },
  "RIGHT": function() { return RIGHT; },
  "ROMAN": function() { return ROMAN; },
  "ROUND": function() { return ROUND; },
  "ROUNDDOWN": function() { return ROUNDDOWN; },
  "ROUNDUP": function() { return ROUNDUP; },
  "ROW": function() { return ROW; },
  "ROWS": function() { return ROWS; },
  "RRI": function() { return RRI; },
  "RSQ": function() { return RSQ; },
  "SEARCH": function() { return SEARCH; },
  "SEC": function() { return SEC; },
  "SECH": function() { return SECH; },
  "SECOND": function() { return SECOND; },
  "SERIESSUM": function() { return SERIESSUM; },
  "SHEET": function() { return SHEET; },
  "SHEETS": function() { return SHEETS; },
  "SIGN": function() { return SIGN; },
  "SIN": function() { return SIN; },
  "SINH": function() { return SINH; },
  "SKEW": function() { return SKEW; },
  "SKEWP": function() { return SKEWP; },
  "SLN": function() { return SLN; },
  "SLOPE": function() { return SLOPE; },
  "SMALL": function() { return SMALL; },
  "SORT": function() { return SORT; },
  "SQRT": function() { return SQRT; },
  "SQRTPI": function() { return SQRTPI; },
  "STANDARDIZE": function() { return STANDARDIZE; },
  "STDEV": function() { return STDEV; },
  "STDEVA": function() { return STDEVA; },
  "STDEVP": function() { return STDEVP; },
  "STDEVPA": function() { return STDEVPA; },
  "STDEVS": function() { return STDEVS; },
  "STEYX": function() { return STEYX; },
  "SUBSTITUTE": function() { return SUBSTITUTE; },
  "SUBTOTAL": function() { return SUBTOTAL; },
  "SUM": function() { return SUM; },
  "SUMIF": function() { return SUMIF; },
  "SUMIFS": function() { return SUMIFS; },
  "SUMPRODUCT": function() { return SUMPRODUCT; },
  "SUMSQ": function() { return SUMSQ; },
  "SUMX2MY2": function() { return SUMX2MY2; },
  "SUMX2PY2": function() { return SUMX2PY2; },
  "SUMXMY2": function() { return SUMXMY2; },
  "SWITCH": function() { return SWITCH; },
  "SYD": function() { return SYD; },
  "T": function() { return T; },
  "TAN": function() { return TAN; },
  "TANH": function() { return TANH; },
  "TBILLEQ": function() { return TBILLEQ; },
  "TBILLPRICE": function() { return TBILLPRICE; },
  "TBILLYIELD": function() { return TBILLYIELD; },
  "TDIST": function() { return TDIST; },
  "TDISTRT": function() { return TDISTRT; },
  "TEXT": function() { return TEXT; },
  "TEXTJOIN": function() { return TEXTJOIN; },
  "TIME": function() { return TIME; },
  "TIMEVALUE": function() { return TIMEVALUE; },
  "TINV": function() { return TINV; },
  "TODAY": function() { return TODAY; },
  "TRANSPOSE": function() { return TRANSPOSE; },
  "TREND": function() { return TREND; },
  "TRIM": function() { return TRIM; },
  "TRIMMEAN": function() { return TRIMMEAN; },
  "TRUE": function() { return TRUE; },
  "TRUNC": function() { return TRUNC; },
  "TTEST": function() { return TTEST; },
  "TYPE": function() { return TYPE; },
  "UNICHAR": function() { return UNICHAR; },
  "UNICODE": function() { return UNICODE; },
  "UNIQUE": function() { return UNIQUE; },
  "UPPER": function() { return UPPER; },
  "VALUE": function() { return VALUE; },
  "VAR": function() { return VAR; },
  "VARA": function() { return VARA; },
  "VARP": function() { return VARP; },
  "VARPA": function() { return VARPA; },
  "VARS": function() { return VARS; },
  "VDB": function() { return VDB; },
  "VLOOKUP": function() { return VLOOKUP; },
  "WEEKDAY": function() { return WEEKDAY; },
  "WEEKNUM": function() { return WEEKNUM; },
  "WEIBULL": function() { return WEIBULL; },
  "WEIBULLDIST": function() { return WEIBULLDIST; },
  "WORKDAY": function() { return WORKDAY; },
  "WORKDAYINTL": function() { return WORKDAYINTL; },
  "XIRR": function() { return XIRR; },
  "XNPV": function() { return XNPV; },
  "XOR": function() { return XOR; },
  "YEAR": function() { return YEAR; },
  "YEARFRAC": function() { return YEARFRAC; },
  "YIELD": function() { return YIELD; },
  "YIELDDISC": function() { return YIELDDISC; },
  "YIELDMAT": function() { return YIELDMAT; },
  "Z": function() { return Z; },
  "ZTEST": function() { return ZTEST; },
  "utils": function() { return utils; }
});

// EXTERNAL MODULE: ./node_modules/jstat/dist/jstat.js
var jstat = __webpack_require__(162);
// EXTERNAL MODULE: ./node_modules/bessel/bessel.js
var bessel = __webpack_require__(765);
;// CONCATENATED MODULE: ./node_modules/@formulajs/formulajs/lib/esm/index.mjs



const nil = new Error('#NULL!');
const div0 = new Error('#DIV/0!');
const value = new Error('#VALUE!');
const ref = new Error('#REF!');
const esm_name = new Error('#NAME?');
const num = new Error('#NUM!');
const na = new Error('#N/A');
const error = new Error('#ERROR!');
const data = new Error('#GETTING_DATA');

var errors = /*#__PURE__*/Object.freeze({
  __proto__: null,
  data: data,
  div0: div0,
  error: error,
  na: na,
  name: esm_name,
  nil: nil,
  num: num,
  ref: ref,
  value: value
});

// Arrays
function argsToArray(args) {
  const result = [];

  arrayEach(args, (value) => {
    result.push(value);
  });

  return result
}

function arrayEach(array, iteratee) {
  let index = -1;
  const length = array.length;

  while (++index < length) {
    if (iteratee(array[index], index, array) === false) {
      break
    }
  }

  return array
}

function arrayValuesToNumbers(arr) {
  let n = arr.length;
  let el;

  while (n--) {
    el = arr[n];

    if (typeof el === 'number') {
      continue
    }

    if (el === true) {
      arr[n] = 1;
      continue
    }

    if (el === false) {
      arr[n] = 0;
      continue
    }

    if (typeof el === 'string') {
      const number = parseNumber(el);

      arr[n] = number instanceof Error ? 0 : number;
    }
  }

  return arr
}

function fillMatrix(matrix, fill_value) {
  if (!matrix) {
    return value
  }

  if (!matrix.every((el) => Array.isArray(el)) || matrix.length === 0) {
    matrix = [[...matrix]];
  }

  matrix.map((arr, i) => {
    arr.map((a, j) => {
      if (!a) {
        matrix[i][j] = 0;
      }
    });
  });

  const longestArrayIndex = matrix.reduce((acc, arr, i) => (arr.length > matrix[acc].length ? i : acc), 0);
  const longestArrayLength = matrix[longestArrayIndex].length;

  return matrix.map((el) => [...el, ...Array(longestArrayLength - el.length).fill(fill_value ? fill_value : 0)])
}

function flatten() {
  let result;

  if (arguments.length === 1) {
    const argument = arguments[0];
    result = isArrayLike(argument) ? argsToArray.apply(null, arguments) : [argument];
  } else {
    result = Array.from(arguments);
  }

  while (!isFlat(result)) {
    result = flattenShallow(result);
  }

  return result
}

function flattenShallow(array) {
  if (!array || !array.reduce) {
    return [array]
  }

  return array.reduce((a, b) => {
    const aIsArray = Array.isArray(a);
    const bIsArray = Array.isArray(b);

    if (aIsArray && bIsArray) {
      return a.concat(b)
    }

    if (aIsArray) {
      a.push(b);

      return a
    }

    if (bIsArray) {
      return [a].concat(b)
    }

    return [a, b]
  })
}

function initial(array, idx) {
  idx = idx || 1;

  if (!array || typeof array.slice !== 'function') {
    return array
  }

  return array.slice(0, array.length - idx)
}

function isArrayLike(a) {
  return a != null && typeof a.length === 'number' && typeof a !== 'string'
}

function isFlat(array) {
  if (!array) {
    return false
  }

  for (let i = 0; i < array.length; ++i) {
    if (Array.isArray(array[i])) {
      return false
    }
  }

  return true
}

function rest(array, idx) {
  idx = idx || 1;

  if (!array || typeof array.slice !== 'function') {
    return array
  }

  return array.slice(idx)
}

function transpose(matrix) {
  if (!matrix) {
    return value
  }

  return matrix[0].map((col, i) => matrix.map((row) => row[i]))
}

// Databases
function findField(database, title) {
  let index = null;

  arrayEach(database, (value, i) => {
    if (value[0] === title) {
      index = i;

      return false
    }
  });

  // Return error if the input field title is incorrect
  if (index == null) {
    return value
  }

  return index
}

// Errors
function anyError() {
  for (let n = 0; n < arguments.length; n++) {
    if (arguments[n] instanceof Error) {
      return arguments[n]
    }
  }

  return undefined
}

function anyIsError() {
  let n = arguments.length;

  while (n--) {
    if (arguments[n] instanceof Error) {
      return true
    }
  }

  return false
}

// Numbers
function cleanFloat(number) {
  const power = 1e14;

  return Math.round(number * power) / power
}

function numbers() {
  const possibleNumbers = flatten.apply(null, arguments);

  return possibleNumbers.filter((el) => typeof el === 'number')
}

function serialNumberToDate(serial) {
  if (serial < 60) {
    serial += 1;
  }

  const utc_days = Math.floor(serial - 25569);
  const utc_value = utc_days * 86400;
  const date_info = new Date(utc_value * 1000);
  const fractional_day = serial - Math.floor(serial) + 0.0000001;

  let total_seconds = Math.floor(86400 * fractional_day);

  const seconds = total_seconds % 60;

  total_seconds -= seconds;

  const hours = Math.floor(total_seconds / (60 * 60));
  const minutes = Math.floor(total_seconds / 60) % 60;
  let days = date_info.getUTCDate();
  let month = date_info.getUTCMonth();

  if (serial >= 60 && serial < 61) {
    days = 29;
    month = 1;
  }

  return new Date(date_info.getUTCFullYear(), month, days, hours, minutes, seconds)
}

// Parsers
function parseBool(bool) {
  if (typeof bool === 'boolean') {
    return bool
  }

  if (bool instanceof Error) {
    return bool
  }

  if (typeof bool === 'number') {
    return bool !== 0
  }

  if (typeof bool === 'string') {
    const up = bool.toUpperCase();

    if (up === 'TRUE') {
      return true
    }

    if (up === 'FALSE') {
      return false
    }
  }

  if (bool instanceof Date && !isNaN(bool)) {
    return true
  }

  return value
}

function parseDate(date) {
  if (!isNaN(date)) {
    if (date instanceof Date) {
      return new Date(date)
    }

    const d = parseFloat(date);

    if (d < 0 || d >= 2958466) {
      return num
    }

    return serialNumberToDate(d)
  }

  if (typeof date === 'string') {
    date = /(\d{4})-(\d\d?)-(\d\d?)$/.test(date) ? new Date(date + 'T00:00:00.000') : new Date(date);

    if (!isNaN(date)) {
      return date
    }
  }

  return value
}

function parseDateArray(arr) {
  let len = arr.length;
  let parsed;

  while (len--) {
    parsed = parseDate(arr[len]);

    if (parsed === value) {
      return parsed
    }

    arr[len] = parsed;
  }

  return arr
}

function parseNumber(string) {
  if (string instanceof Error) {
    return string
  }

  if (string === undefined || string === null) {
    return 0
  }

  if (typeof string === 'boolean') {
    string = +string;
  }

  if (!isNaN(string) && string !== '') {
    return parseFloat(string)
  }

  return value
}

function parseNumberArray(arr) {
  let len;

  if (!arr || (len = arr.length) === 0) {
    return value
  }

  let parsed;

  while (len--) {
    if (arr[len] instanceof Error) {
      return arr[len]
    }

    parsed = parseNumber(arr[len]);

    if (parsed instanceof Error) {
      return parsed
    }

    arr[len] = parsed;
  }

  return arr
}

function parseString(string) {
  if (string instanceof Error) {
    return string
  }

  if (string === undefined || string === null) {
    return ''
  }

  return string.toString()
}

// Strings
function anyIsString() {
  let n = arguments.length;

  while (n--) {
    if (typeof arguments[n] === 'string') {
      return true
    }
  }

  return false
}

// Misc
function isDefined(arg) {
  return arg !== undefined && arg !== null
}

const defaultOperator = '=';
const validSymbols = ['>', '>=', '<', '<=', '=', '<>'];
const _TOKEN_TYPE_OPERATOR = 'operator';
const _TOKEN_TYPE_LITERAL = 'literal';
const SUPPORTED_TOKENS = [_TOKEN_TYPE_OPERATOR, _TOKEN_TYPE_LITERAL];

const TOKEN_TYPE_OPERATOR = _TOKEN_TYPE_OPERATOR;
const TOKEN_TYPE_LITERAL = _TOKEN_TYPE_LITERAL;

/**
 * Create token which describe passed symbol/value.
 *
 * @param {String} value Value/Symbol to describe.
 * @param {String} type Type of the token 'operator' or 'literal'.
 * @return {Object}
 */
function createToken(value, type) {
  if (SUPPORTED_TOKENS.indexOf(type) === -1) {
    throw new Error('Unsupported token type: ' + type)
  }

  return {
    value: value,
    type: type
  }
}

/**
 * Tries to cast numeric values to their type passed as a string.
 *
 * @param {*} value
 * @return {*}
 */
function castValueToCorrectType(value) {
  if (typeof value !== 'string') {
    return value
  }

  if (/^\d+(\.\d+)?$/.test(value)) {
    value = value.indexOf('.') === -1 ? parseInt(value, 10) : parseFloat(value);
  }

  return value
}

/**
 * Generate stream of tokens from passed expression.
 *
 * @param {String} expression
 * @return {String[]}
 */
function tokenizeExpression(expression) {
  const expressionLength = expression.length;
  const tokens = [];
  let cursorIndex = 0;
  let processedValue = '';
  let processedSymbol = '';

  while (cursorIndex < expressionLength) {
    const char = expression.charAt(cursorIndex);

    switch (char) {
      case '>':
      case '<':
      case '=':
        processedSymbol = processedSymbol + char;

        if (processedValue.length > 0) {
          tokens.push(processedValue);
          processedValue = '';
        }

        break
      default:
        if (processedSymbol.length > 0) {
          tokens.push(processedSymbol);
          processedSymbol = '';
        }

        processedValue = processedValue + char;
        break
    }

    cursorIndex++;
  }

  if (processedValue.length > 0) {
    tokens.push(processedValue);
  }

  if (processedSymbol.length > 0) {
    tokens.push(processedSymbol);
  }

  return tokens
}

/**
 * Analyze and convert tokens to an object which describes their meaning.
 *
 * @param {String[]} tokens
 * @return {Object[]}
 */
function analyzeTokens(tokens) {
  let literalValue = '';
  const analyzedTokens = [];

  for (let i = 0; i < tokens.length; i++) {
    const token = tokens[i];

    if (i === 0 && validSymbols.indexOf(token) >= 0) {
      analyzedTokens.push(createToken(token, TOKEN_TYPE_OPERATOR));
    } else {
      literalValue += token;
    }
  }

  if (literalValue.length > 0) {
    analyzedTokens.push(createToken(castValueToCorrectType(literalValue), TOKEN_TYPE_LITERAL));
  }

  if (analyzedTokens.length > 0 && analyzedTokens[0].type !== TOKEN_TYPE_OPERATOR) {
    analyzedTokens.unshift(createToken(defaultOperator, TOKEN_TYPE_OPERATOR));
  }

  return analyzedTokens
}

/**
 * Compute/Evaluate an expression passed as an array of tokens.
 *
 * @param {Object[]} tokens
 * @return {Boolean}
 */
function computeExpression(tokens) {
  const values = [];
  let operator;

  for (let i = 0; i < tokens.length; i++) {
    const token = tokens[i];

    switch (token.type) {
      case TOKEN_TYPE_OPERATOR:
        operator = token.value;
        break
      case TOKEN_TYPE_LITERAL:
        values.push(token.value);
        break
    }
  }

  return evaluate(values, operator)
}

/**
 * Evaluate values based on passed math operator.
 *
 * @param {*} values
 * @param {String} operator
 * @return {Boolean}
 */
function evaluate(values, operator) {
  let result = false;

  switch (operator) {
    case '>':
      result = values[0] > values[1];
      break
    case '>=':
      result = values[0] >= values[1];
      break
    case '<':
      result = values[0] < values[1];
      break
    case '<=':
      result = values[0] <= values[1];
      break
    case '=':
      result = values[0] == values[1];
      break
    case '<>':
      result = values[0] != values[1];
      break
  }

  return result
}

function parse(expression) {
  return analyzeTokens(tokenizeExpression(expression))
}

const compute = computeExpression;

// TODO
/**
 * -- Not implemented --
 *
 * Returns information about the formatting, location, or contents of a value.
 *
 * Category: Information
 *
 * @returns
 */
function CELL() {
  throw new Error('CELL is not implemented')
}

const ERROR = {};

ERROR.TYPE = (error_val) => {
  switch (error_val) {
    case nil:
      return 1
    case div0:
      return 2
    case value:
      return 3
    case ref:
      return 4
    case esm_name:
      return 5
    case num:
      return 6
    case na:
      return 7
    case data:
      return 8
  }

  return na
};

// TODO
/**
 * -- Not implemented --
 *
 * Returns information about the current operating environment.
 *
 * Category: Information
 *
 * @returns
 */
function INFO() {
  throw new Error('INFO is not implemented')
}

/**
 * Returns TRUE if the value is blank.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISBLANK(value) {
  return value === null
}

/**
 * Returns TRUE if the value is any error value except #N/A.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISERR(value$1) {
  return (
    [value, ref, div0, num, esm_name, nil].indexOf(value$1) >= 0 ||
    (typeof value$1 === 'number' && (isNaN(value$1) || !isFinite(value$1)))
  )
}

/**
 * Returns TRUE if the value is any error value.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISERROR(value) {
  return ISERR(value) || value === na
}

/**
 * Returns TRUE if the number is even.
 *
 * Category: Information
 *
 * @param {*} number The value to test. If number is not an integer, it is truncated.
 * @returns
 */
function ISEVEN(number) {
  return !(Math.floor(Math.abs(number)) & 1)
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns TRUE if there is a reference to a value that contains a formula.
 *
 * Category: Information
 *
 * @param {*} reference Reference is a reference to the value you want to test. Reference can be a value reference, a formula, or a name that refers to a value.
 * @returns
 */
function ISFORMULA() {
  throw new Error('ISFORMULA is not implemented')
}

/**
 * Returns TRUE if the value is a logical value.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISLOGICAL(value) {
  return value === true || value === false
}

/**
 * Returns TRUE if the value is the #N/A error value.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISNA(value) {
  return value === na
}

/**
 * Returns TRUE if the value is not text.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISNONTEXT(value) {
  return typeof value !== 'string'
}

/**
 * Returns TRUE if the value is a number.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISNUMBER(value) {
  return typeof value === 'number' && !isNaN(value) && isFinite(value)
}

/**
 * Returns TRUE if the number is odd.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISODD(value) {
  return !!(Math.floor(Math.abs(value)) & 1)
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns TRUE if the value is a reference.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISREF() {
  throw new Error('ISREF is not implemented')
}

/**
 * Returns TRUE if the value is text.
 *
 * Category: Information
 *
 * @param {*} value The value that you want tested. The value argument can be a blank (empty value), error, logical value, text, number, or reference value, or a name referring to any of these.
 * @returns
 */
function ISTEXT(value) {
  return typeof value === 'string'
}

/**
 * Returns a value converted to a number.
 *
 * Category: Information
 *
 * @param {*} value The value you want converted. N converts values listed in the following table.
 * @returns
 */
function N(value) {
  if (ISNUMBER(value)) {
    return value
  }

  if (value instanceof Date) {
    return value.getTime()
  }

  if (value === true) {
    return 1
  }

  if (value === false) {
    return 0
  }

  if (ISERROR(value)) {
    return value
  }

  return 0
}

/**
 * Returns the error value #N/A.
 *
 * Category: Information
 *
 * @returns
 */
function NA() {
  return na
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the sheet number of the referenced sheet.
 *
 * Category: Information
 *
 * @param {*} value Optional. Value is the name of a sheet or a reference for which you want the sheet number. If value is omitted, SHEET returns the number of the sheet that contains the function.
 * @returns
 */
function SHEET() {
  throw new Error('SHEET is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the number of sheets in a reference.
 *
 * Category: Information
 *
 * @param {*} reference Optional. Reference is a reference for which you want to know the number of sheets it contains. If Reference is omitted, SHEETS returns the number of sheets in the workbook that contains the function.
 * @returns
 */
function SHEETS() {
  throw new Error('SHEETS is not implemented')
}

/**
 * Returns a number indicating the data type of a value.
 *
 * Category: Information
 *
 * @param {*} value Can be any Microsoft Excel value, such as a number, text, logical value, and so on.
 * @returns
 */
function TYPE(value) {
  if (ISNUMBER(value)) {
    return 1
  }

  if (ISTEXT(value)) {
    return 2
  }

  if (ISLOGICAL(value)) {
    return 4
  }

  if (ISERROR(value)) {
    return 16
  }

  if (Array.isArray(value)) {
    return 64
  }
}

/**
 * Chooses a value from a list of values.
 *
 * Category: Lookup and reference
 *
 * @param {*} index_num Specifies which value argument is selected. Index_num must be a number between 1 and 254, or a formula or reference to a value containing a number between 1 and 254. If index_num is 1, CHOOSE returns value1; if it is 2, CHOOSE returns value2; and so on. If index_num is less than 1 or greater than the number of the last value in the list, CHOOSE returns the #VALUE! error value. If index_num is a fraction, it is truncated to the lowest integer before being used.
 - If index_num is 1, CHOOSE returns value1; if it is 2, CHOOSE returns value2; and so on.
 - If index_num is less than 1 or greater than the number of the last value in the list, CHOOSE returns the #VALUE! error value.
 - If index_num is a fraction, it is truncated to the lowest integer before being used.
 * @param {*} args value1, value2, ... Value 1 is required, subsequent values are optional. 1 to 254 value arguments from which CHOOSE selects a value or an action to perform based on index_num. The arguments can be numbers, value references, defined names, formulas, functions, or text.
 * @returns
 */
function CHOOSE() {
  if (arguments.length < 2) {
    return na
  }

  const index = arguments[0];

  if (index < 1 || index > 254) {
    return value
  }

  if (arguments.length < index + 1) {
    return value
  }

  return arguments[index]
}

/**
 * Returns the column number of a reference.
 *
 * Category: Lookup and reference
 *
 * @param {*} reference the value or range of values for which you want to return the column number.
 * @param {*} index
 * @returns
 */
function COLUMN(reference, index) {
  if (arguments.length !== 2) {
    return na
  }

  if (index < 0) {
    return num
  }

  if (!(reference instanceof Array) || typeof index !== 'number') {
    return value
  }

  if (reference.length === 0) {
    return undefined
  }

  return jstat.col(reference, index)
}

/**
 * Returns the number of columns in a reference.
 *
 * Category: Lookup and reference
 *
 * @param {*} array An array or array formula, or a reference to a range of values for which you want the number of columns.
 * @returns
 */
function COLUMNS(array) {
  if (arguments.length !== 1) {
    return na
  }

  if (!(array instanceof Array)) {
    return value
  }

  if (array.length === 0) {
    return 0
  }

  return jstat.cols(array)
}

/**
 * Looks in the top row of an array and returns the value of the indicated value.
 *
 * Category: Lookup and reference
 *
 * @param {*} lookup_value The value to be found in the first row of the table. Lookup_value can be a value, a reference, or a text string.
 * @param {*} table_array A table of information in which data is looked up. Use a reference to a range or a range name.
 * @param {*} row_index_num The row number in table_array from which the matching value will be returned. A row_index_num of 1 returns the first row value in table_array, a row_index_num of 2 returns the second row value in table_array, and so on. If row_index_num is less than 1, HLOOKUP returns the #VALUE! error value; if row_index_num is greater than the number of rows on table_array, HLOOKUP returns the #REF! error value.
 * @param {*} range_lookup Optional. A logical value that specifies whether you want HLOOKUP to find an exact match or an approximate match. If TRUE or omitted, an approximate match is returned. In other words, if an exact match is not found, the next largest value that is less than lookup_value is returned. If FALSE, HLOOKUP will find an exact match. If one is not found, the error value #N/A is returned.
 * @returns
 */
function HLOOKUP(lookup_value, table_array, row_index_num, range_lookup) {
  return VLOOKUP(lookup_value, transpose(table_array), row_index_num, range_lookup)
}

/**
 * Uses an index to choose a value from a reference or array.
 *
 * Category: Lookup and reference
 *
 * @param {*} array A range of values or an array constant.
 - If array contains only one row or column, the corresponding row_num or column_num argument is optional.
 - If array has more than one row and more than one column, and only row_num or column_num is used, INDEX returns an array of the entire row or column in array.
 * @param {*} row_num Required, unless column_num is present. Selects the row in array from which to return a value. If row_num is omitted, column_num is required.
 * @param {*} column_num Optional. Selects the column in array from which to return a value. If column_num is omitted, row_num is required.
 * @returns
 */
function INDEX(array, row_num, column_num) {
  const someError = anyError(array, row_num, column_num);

  if (someError) {
    return someError
  }

  if (!Array.isArray(array)) {
    return value
  }

  const isOneDimensionRange = array.length > 0 && !Array.isArray(array[0]);

  if (isOneDimensionRange && !column_num) {
    column_num = row_num;
    row_num = 1;
  } else {
    column_num = column_num || 1;
    row_num = row_num || 1;
  }

  if (column_num < 0 || row_num < 0) {
    return value
  }

  if (isOneDimensionRange && row_num === 1 && column_num <= array.length) {
    return array[column_num - 1]
  } else if (row_num <= array.length && column_num <= array[row_num - 1].length) {
    return array[row_num - 1][column_num - 1]
  }

  return ref
}

/**
 * Looks up values in a vector or array.
 *
 * Category: Lookup and reference
 *
 * @param {*} lookup_value A value that LOOKUP searches for in an array. The lookup_value argument can be a number, text, a logical value, or a name or reference that refers to a value.
 - If LOOKUP can't find the value of lookup_value, it uses the largest value in the array that is less than or equal to lookup_value.
 - If the value of lookup_value is smaller than the smallest value in the first row or column (depending on the array dimensions), LOOKUP returns the #N/A error value.
 * @param {*} array A range of values that contains text, numbers, or logical values that you want to compare with lookup_value. The array form of LOOKUP is very similar to the HLOOKUP and VLOOKUP functions. The difference is that HLOOKUP searches for the value of lookup_value in the first row, VLOOKUP searches in the first column, and LOOKUP searches according to the dimensions of array.
* @param {*} result_array Optional. A range that contains only one row or column. The result_array argument must be the same size as lookup_value. It has to be the same size.
 * @returns
 */
function LOOKUP(lookup_value, array, result_array) {
  array = flatten(array);
  result_array = result_array ? flatten(result_array) : array;

  const isNumberLookup = typeof lookup_value === 'number';
  let result = na;

  for (let i = 0; i < array.length; i++) {
    if (array[i] === lookup_value) {
      return result_array[i]
    } else if (
      (isNumberLookup && array[i] <= lookup_value) ||
      (typeof array[i] === 'string' && array[i].localeCompare(lookup_value) < 0)
    ) {
      result = result_array[i];
    } else if (isNumberLookup && array[i] > lookup_value) {
      return result
    }
  }

  return result
}

/**
 * Looks up values in a reference or array.
 *
 * Category: Lookup and reference
 *
 * @param {*} lookup_value The value that you want to match in lookup_array. For example, when you look up someone's number in a telephone book, you are using the person's name as the lookup value, but the telephone number is the value you want.The lookup_value argument can be a value (number, text, or logical value) or a value reference to a number, text, or logical value.
 * @param {*} lookup_array The range of values being searched.
 * @param {*} match_type Optional. The number -1, 0, or 1. The match_type argument specifies how Excel matches lookup_value with values in lookup_array. The default value for this argument is 1.
 * @returns
 */
function MATCH(lookup_value, lookup_array, match_type) {
  if (!lookup_value || !lookup_array) {
    return na
  }

  if (arguments.length === 2) {
    match_type = 1;
  }

  lookup_array = flatten(lookup_array);

  if (!(lookup_array instanceof Array)) {
    return na
  }

  if (match_type !== -1 && match_type !== 0 && match_type !== 1) {
    return na
  }

  let index;
  let indexValue;

  for (let idx = 0; idx < lookup_array.length; idx++) {
    if (match_type === 1) {
      if (lookup_array[idx] === lookup_value) {
        return idx + 1
      } else if (lookup_array[idx] < lookup_value) {
        if (!indexValue) {
          index = idx + 1;
          indexValue = lookup_array[idx];
        } else if (lookup_array[idx] > indexValue) {
          index = idx + 1;
          indexValue = lookup_array[idx];
        }
      }
    } else if (match_type === 0) {
      if (typeof lookup_value === 'string' && typeof lookup_array[idx] === 'string') {
        const lookupValueStr = lookup_value.toLowerCase().replace(/\?/g, '.').replace(/\*/g, '.*').replace(/~/g, '\\');
        const regex = new RegExp('^' + lookupValueStr + '$');

        if (regex.test(lookup_array[idx].toLowerCase())) {
          return idx + 1
        }
      } else {
        if (lookup_array[idx] === lookup_value) {
          return idx + 1
        }
      }
    } else if (match_type === -1) {
      if (lookup_array[idx] === lookup_value) {
        return idx + 1
      } else if (lookup_array[idx] > lookup_value) {
        if (!indexValue) {
          index = idx + 1;
          indexValue = lookup_array[idx];
        } else if (lookup_array[idx] < indexValue) {
          index = idx + 1;
          indexValue = lookup_array[idx];
        }
      }
    }
  }

  return index || na
}

/**
 * Returns the number of rows in a reference.
 *
 * Category: Lookup and reference
 *
 * @param {*} array An array, an array formula, or a reference to a range of values for which you want the number of rows.
 * @returns
 */
function ROWS(array) {
  if (arguments.length !== 1) {
    return na
  }

  if (!(array instanceof Array)) {
    return value
  }

  if (array.length === 0) {
    return 0
  }

  return jstat.rows(array)
}
/**
 * Returns a sorted array of the elements in an array. The returned array is the same shape as the provided array argument.
 *
 * Category: Lookup and reference
 *
 * @param {*} array Array to sort
 * @param {*} sort_index Optional. A number indicating the row or column to sort by
 * @param {*} sort_order Optional. A number indicating the desired sort order; 1 for ascending order (default), -1 for descending order
 * @param {*} by_col Optional. A logical value indicating the desired sort direction; FALSE to sort by row (default), TRUE to sort by column
 * @returns
 */
function SORT(array, sort_index = 1, sort_order = 1, by_col = false) {
  if (!array || !Array.isArray(array)) {
    return na
  }

  if (array.length === 0) {
    return 0
  }

  sort_index = parseNumber(sort_index);
  if (!sort_index || sort_index < 1) {
    return value
  }

  sort_order = parseNumber(sort_order);
  if (sort_order !== 1 && sort_order !== -1) {
    return value
  }

  by_col = parseBool(by_col);
  if (typeof by_col !== 'boolean') {
    return esm_name
  }

  const sortArray = (arr) =>
    arr.sort((a, b) => {
      a = parseString(a[sort_index - 1]);
      b = parseString(b[sort_index - 1]);

      return sort_order === 1 ? (a < b ? sort_order * -1 : sort_order) : a > b ? sort_order : sort_order * -1
    });

  const matrix = fillMatrix(array);
  const result = by_col ? transpose(matrix) : matrix;

  return sort_index >= 1 && sort_index <= result[0].length
    ? by_col
      ? transpose(sortArray(result))
      : sortArray(result)
    : value
}

/**
 * Returns the transpose of an array.
 *
 * Category: Lookup and reference
 *
 * @param {*} array An array or range of values on a worksheet that you want to transpose. The transpose of an array is created by using the first row of the array as the first column of the new array, the second row of the array as the second column of the new array, and so on. If you're not sure of how to enter an array formula, see Create an array formula.
 * @returns
 */
function TRANSPOSE(array) {
  if (!array) {
    return na
  }

  const matrix = fillMatrix(array);

  return transpose(matrix)
}

/**
 * Returns a list of unique values in a list or range.
 *
 * Category: Lookup and reference
 *
 * @returns
 */
function UNIQUE() {
  const result = [];

  for (let i = 0; i < arguments.length; ++i) {
    let hasElement = false;
    const element = arguments[i];

    // Check if we've already seen this element.

    for (let j = 0; j < result.length; ++j) {
      hasElement = result[j] === element;

      if (hasElement) {
        break
      }
    }

    // If we did not find it, add it to the result.
    if (!hasElement) {
      result.push(element);
    }
  }

  return result
}

/**
 * Looks in the first column of an array and moves across the row to return the value of a value.
 *
 * Category: Lookup and reference
 *
 * @param {*} lookup_value The value to be found in the first row of the table. Lookup_value can be a value, a reference, or a text string.
 * @param {*} table_array A table of information in which data is looked up. Use a reference to a range or a range name.
 * @param {*} col_index_num The row number in table_array from which the matching value will be returned. A row_index_num of 1 returns the first row value in table_array, a row_index_num of 2 returns the second row value in table_array, and so on. If row_index_num is less than 1, HLOOKUP returns the #VALUE! error value; if row_index_num is greater than the number of rows on table_array, HLOOKUP returns the #REF! error value.
 * @param {*} range_lookup Optional. A logical value that specifies whether you want HLOOKUP to find an exact match or an approximate match. If TRUE or omitted, an approximate match is returned. In other words, if an exact match is not found, the next largest value that is less than lookup_value is returned. If FALSE, HLOOKUP will find an exact match. If one is not found, the error value #N/A is returned.
 * @returns
 */
function VLOOKUP(lookup_value, table_array, col_index_num, range_lookup) {
  if (!table_array || !col_index_num) {
    return na
  }

  range_lookup = !(range_lookup === 0 || range_lookup === false);
  let result = na;
  const isNumberLookup = typeof lookup_value === 'number';
  let exactMatchOnly = false;

  for (let i = 0; i < table_array.length; i++) {
    const row = table_array[i];

    if (row[0] === lookup_value) {
      result = col_index_num < row.length + 1 ? row[col_index_num - 1] : ref;
      break
    } else if (
      !exactMatchOnly &&
      ((isNumberLookup && range_lookup && row[0] <= lookup_value) ||
        (range_lookup && typeof row[0] === 'string' && row[0].localeCompare(lookup_value) < 0))
    ) {
      result = col_index_num < row.length + 1 ? row[col_index_num - 1] : ref;
    }

    if (isNumberLookup && row[0] > lookup_value) {
      exactMatchOnly = true;
    }
  }

  return result
}

const SQRT2PI = 2.5066282746310002;

/**
 * Returns the average of the absolute deviations of data points from their mean.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 arguments for which you want the average of the absolute deviations. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
function AVEDEV() {
  const flatArguments = flatten(arguments);
  const flatArgumentsDefined = flatArguments.filter(isDefined);

  if (flatArgumentsDefined.length === 0) {
    return num
  }

  const range = parseNumberArray(flatArgumentsDefined);

  if (range instanceof Error) {
    return range
  }

  return jstat.sum(jstat(range).subtract(jstat.mean(range)).abs()[0]) / range.length
}

/**
 * Returns the average of its arguments.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ...Numbers, value references or ranges for which you want the average.
 * @returns
 */
function AVERAGE() {
  const flatArguments = flatten(arguments);
  const flatArgumentsDefined = flatArguments.filter(isDefined);

  if (flatArgumentsDefined.length === 0) {
    return div0
  }

  const someError = anyError.apply(undefined, flatArgumentsDefined);

  if (someError) {
    return someError
  }

  const range = numbers(flatArgumentsDefined);
  const n = range.length;
  let sum = 0;
  let count = 0;
  let result;

  for (let i = 0; i < n; i++) {
    sum += range[i];
    count += 1;
  }

  result = sum / count;

  if (isNaN(result)) {
    result = num;
  }

  return result
}

/**
 * Returns the average of its arguments, including numbers, text, and logical values.
 *
 * Category: Statistical
 *
 * @param {*} args value1, value2, ... Value1 is required, subsequent values are optional. 1 to 255 values, ranges of values, or values for which you want the average.
 * @returns
 */
function AVERAGEA() {
  const flatArguments = flatten(arguments);
  const flatArgumentsDefined = flatArguments.filter(isDefined);

  if (flatArgumentsDefined.length === 0) {
    return div0
  }

  const someError = anyError.apply(undefined, flatArgumentsDefined);

  if (someError) {
    return someError
  }

  const range = flatArgumentsDefined;
  const n = range.length;
  let sum = 0;
  let count = 0;
  let result;

  for (let i = 0; i < n; i++) {
    const el = range[i];

    if (typeof el === 'number') {
      sum += el;
    }

    if (el === true) {
      sum++;
    }

    if (el !== null) {
      count++;
    }
  }

  result = sum / count;

  if (isNaN(result)) {
    result = num;
  }

  return result
}

/**
 * Returns the average (arithmetic mean) of all the values in a range that meet a given criteria.
 *
 * Category: Statistical
 *
 * @param {*} range One or more values to average, including numbers or names, arrays, or references that contain numbers.
 * @param {*} criteria The criteria in the form of a number, expression, value reference, or text that defines which values are averaged.
 * @param {*} average_range Optional. The actual set of values to average. If omitted, range is used.
 * @returns
 */
function AVERAGEIF(range, criteria, average_range) {
  if (arguments.length <= 1) {
    return na
  }

  average_range = average_range || range;
  const flatAverageRange = flatten(average_range);
  const flatAverageRangeDefined = flatAverageRange.filter(isDefined);
  average_range = parseNumberArray(flatAverageRangeDefined);

  range = flatten(range);

  if (average_range instanceof Error) {
    return average_range
  }

  let average_count = 0;
  let result = 0;
  const isWildcard = criteria === void 0 || criteria === '*';
  const tokenizedCriteria = isWildcard ? null : parse(criteria + '');

  for (let i = 0; i < range.length; i++) {
    const value = range[i];

    if (isWildcard) {
      result += average_range[i];
      average_count++;
    } else {
      const tokens = [createToken(value, TOKEN_TYPE_LITERAL)].concat(tokenizedCriteria);

      if (compute(tokens)) {
        result += average_range[i];
        average_count++;
      }
    }
  }

  return result / average_count
}

/**
 * Returns the average (arithmetic mean) of all values that meet multiple criteria.
 *
 * Category: Statistical
 *
 * @param {*} args One or more values to average, including numbers or names, arrays, or references that contain numbers.
 * @returns
 */
function AVERAGEIFS() {
  // Does not work with multi dimensional ranges yet!
  // http://office.microsoft.com/en-001/excel-help/averageifs-function-HA010047493.aspx
  const args = argsToArray(arguments);
  const criteriaLength = (args.length - 1) / 2;
  const range = flatten(args[0]);
  let count = 0;
  let result = 0;

  for (let i = 0; i < range.length; i++) {
    let isMeetCondition = false;

    for (let j = 0; j < criteriaLength; j++) {
      const value = args[2 * j + 1][i];
      const criteria = args[2 * j + 2];
      const isWildcard = criteria === void 0 || criteria === '*';
      let computedResult = false;

      if (isWildcard) {
        computedResult = true;
      } else {
        const tokenizedCriteria = parse(criteria + '');
        const tokens = [createToken(value, TOKEN_TYPE_LITERAL)].concat(tokenizedCriteria);

        computedResult = compute(tokens);
      }

      // Criterias are calculated as AND so any `false` breakes the loop as unmeet condition
      if (!computedResult) {
        isMeetCondition = false;
        break
      }

      isMeetCondition = true;
    }

    if (isMeetCondition) {
      result += range[i];
      count++;
    }
  }

  const average = result / count;

  return isNaN(average) ? 0 : average
}

const BETA = {};

/**
 * Returns the beta cumulative distribution function.
 *
 * Category: Statistical
 *
 * @param {*} x The value between A and B at which to evaluate the function
 * @param {*} alpha A parameter of the distribution.
 * @param {*} beta A parameter of the distribution.
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, BETA.DIST returns the cumulative distribution function; if FALSE, it returns the probability density function.
 * @param {*} a Optional. A lower bound to the interval of x.
 * @param {*} b Optional. An upper bound to the interval of x.
 * @returns
 */
BETA.DIST = function (x, alpha, beta, cumulative, a, b) {
  if (arguments.length < 4) {
    return value
  }

  a = a === undefined ? 0 : a;
  b = b === undefined ? 1 : b;

  x = parseNumber(x);
  alpha = parseNumber(alpha);
  beta = parseNumber(beta);
  a = parseNumber(a);
  b = parseNumber(b);

  if (anyIsError(x, alpha, beta, a, b)) {
    return value
  }

  x = (x - a) / (b - a);

  return cumulative ? jstat.beta.cdf(x, alpha, beta) : jstat.beta.pdf(x, alpha, beta)
};

/**
 * Returns the inverse of the cumulative distribution function for a specified beta distribution.
 *
 * Category: Statistical
 *
 * @param {*} probability A probability associated with the beta distribution.
 * @param {*} alpha A parameter of the distribution.
 * @param {*} beta A parameter the distribution.
 * @param {*} a Optional. A lower bound to the interval of x.
 * @param {*} b Optional. An upper bound to the interval of x.
 * @returns
 */
BETA.INV = (probability, alpha, beta, a, b) => {
  a = a === undefined ? 0 : a;
  b = b === undefined ? 1 : b;
  probability = parseNumber(probability);
  alpha = parseNumber(alpha);
  beta = parseNumber(beta);
  a = parseNumber(a);
  b = parseNumber(b);

  if (anyIsError(probability, alpha, beta, a, b)) {
    return value
  }

  return jstat.beta.inv(probability, alpha, beta) * (b - a) + a
};

const BINOM = {};

/**
 * Returns the individual term binomial distribution probability.
 *
 * Category: Statistical
 *
 * @param {*} number_s The number of successes in trials.
 * @param {*} trials The number of independent trials.
 * @param {*} probability_s The probability of success on each trial.
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, then BINOM.DIST returns the cumulative distribution function, which is the probability that there are at most number_s successes; if FALSE, it returns the probability mass function, which is the probability that there are number_s successes.
 * @returns
 */
BINOM.DIST = (number_s, trials, probability_s, cumulative) => {
  number_s = parseNumber(number_s);
  trials = parseNumber(trials);
  probability_s = parseNumber(probability_s);
  cumulative = parseNumber(cumulative);

  if (anyIsError(number_s, trials, probability_s, cumulative)) {
    return value
  }

  return cumulative
    ? jstat.binomial.cdf(number_s, trials, probability_s)
    : jstat.binomial.pdf(number_s, trials, probability_s)
};

/**
 * Returns the probability of a trial result using a binomial distribution.
 *
 * Category: Statistical
 *
 * @param {*} trials The number of independent trials. Must be greater than or equal to 0.
 * @param {*} probability_s The probability of success in each trial. Must be greater than or equal to 0 and less than or equal to 1.
 * @param {*} number_s The number of successes in trials. Must be greater than or equal to 0 and less than or equal to Trials.
 * @param {*} number_s2 Optional. If provided, returns the probability that the number of successful trials will fall between Number_s and number_s2. Must be greater than or equal to Number_s and less than or equal to Trials.
 * @returns
 */
BINOM.DIST.RANGE = (trials, probability_s, number_s, number_s2) => {
  number_s2 = number_s2 === undefined ? number_s : number_s2;

  trials = parseNumber(trials);
  probability_s = parseNumber(probability_s);
  number_s = parseNumber(number_s);
  number_s2 = parseNumber(number_s2);

  if (anyIsError(trials, probability_s, number_s, number_s2)) {
    return value
  }

  let result = 0;

  for (let i = number_s; i <= number_s2; i++) {
    result += COMBIN(trials, i) * Math.pow(probability_s, i) * Math.pow(1 - probability_s, trials - i);
  }

  return result
};

/**
 * Returns the smallest value for which the cumulative binomial distribution is less than or equal to a criterion value.
 *
 * Category: Statistical
 *
 * @param {*} trials The number of Bernoulli trials.
 * @param {*} probability_s The probability of a success on each trial.
 * @param {*} alpha The criterion value.
 * @returns
 */
BINOM.INV = (trials, probability_s, alpha) => {
  trials = parseNumber(trials);
  probability_s = parseNumber(probability_s);
  alpha = parseNumber(alpha);

  if (anyIsError(trials, probability_s, alpha)) {
    return value
  }

  let x = 0;

  while (x <= trials) {
    if (jstat.binomial.cdf(x, trials, probability_s) >= alpha) {
      return x
    }

    x++;
  }
};

const CHISQ = {};

/**
 * Returns the cumulative beta probability density function.
 *
 * Category: Statistical
 *
 * @param {*} x The value at which you want to evaluate the distribution.
 * @param {*} deg_freedom The number of degrees of freedom.
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, CHISQ.DIST returns the cumulative distribution function; if FALSE, it returns the probability density function.
 * @returns
 */
CHISQ.DIST = (x, deg_freedom, cumulative) => {
  x = parseNumber(x);
  deg_freedom = parseNumber(deg_freedom);

  if (anyIsError(x, deg_freedom)) {
    return value
  }

  return cumulative ? jstat.chisquare.cdf(x, deg_freedom) : jstat.chisquare.pdf(x, deg_freedom)
};

/**
 * Returns the one-tailed probability of the chi-squared distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The value at which you want to evaluate the distribution.
 * @param {*} deg_freedom The number of degrees of freedom.
 * @returns
 */
CHISQ.DIST.RT = (x, deg_freedom) => {
  if (!x | !deg_freedom) {
    return na
  }

  if (x < 1 || deg_freedom > Math.pow(10, 10)) {
    return num
  }

  if (typeof x !== 'number' || typeof deg_freedom !== 'number') {
    return value
  }

  return 1 - jstat.chisquare.cdf(x, deg_freedom)
};

/**
 * Returns the cumulative beta probability density function.
 *
 * Category: Statistical
 *
 * @param {*} probability A probability associated with the chi-squared distribution.
 * @param {*} deg_freedom The number of degrees of freedom.
 * @returns
 */
CHISQ.INV = (probability, deg_freedom) => {
  probability = parseNumber(probability);
  deg_freedom = parseNumber(deg_freedom);

  if (anyIsError(probability, deg_freedom)) {
    return value
  }

  return jstat.chisquare.inv(probability, deg_freedom)
};

/**
 * Returns the inverse of the one-tailed probability of the chi-squared distribution.
 *
 * Category: Statistical
 *
 * @param {*} probability A probability associated with the chi-squared distribution.
 * @param {*} deg_freedom The number of degrees of freedom.
 * @returns
 */
CHISQ.INV.RT = (probability, deg_freedom) => {
  if (!probability | !deg_freedom) {
    return na
  }

  if (probability < 0 || probability > 1 || deg_freedom < 1 || deg_freedom > Math.pow(10, 10)) {
    return num
  }

  if (typeof probability !== 'number' || typeof deg_freedom !== 'number') {
    return value
  }

  return jstat.chisquare.inv(1.0 - probability, deg_freedom)
};

/**
 * Returns the test for independence.
 *
 * Category: Statistical
 *
 * @param {*} actual_range The range of data that contains observations to test against expected values.
 * @param {*} expected_range The range of data that contains the ratio of the product of row totals and column totals to the grand total.
 * @returns
 */
CHISQ.TEST = function (actual_range, expected_range) {
  if (arguments.length !== 2) {
    return na
  }

  if (!(actual_range instanceof Array) || !(expected_range instanceof Array)) {
    return value
  }

  if (actual_range.length !== expected_range.length) {
    return value
  }

  if (actual_range[0] && expected_range[0] && actual_range[0].length !== expected_range[0].length) {
    return value
  }

  const row = actual_range.length;
  let tmp, i, j;

  // Convert single-dimension array into two-dimension array

  for (i = 0; i < row; i++) {
    if (!(actual_range[i] instanceof Array)) {
      tmp = actual_range[i];
      actual_range[i] = [];
      actual_range[i].push(tmp);
    }

    if (!(expected_range[i] instanceof Array)) {
      tmp = expected_range[i];
      expected_range[i] = [];
      expected_range[i].push(tmp);
    }
  }

  const col = actual_range[0].length;
  const dof = col === 1 ? row - 1 : (row - 1) * (col - 1);
  let xsqr = 0;
  const Pi = Math.PI;

  for (i = 0; i < row; i++) {
    for (j = 0; j < col; j++) {
      xsqr += Math.pow(actual_range[i][j] - expected_range[i][j], 2) / expected_range[i][j];
    }
  }

  // Get independency by X square and its degree of freedom
  function ChiSq(xsqr, dof) {
    let p = Math.exp(-0.5 * xsqr);

    if (dof % 2 === 1) {
      p = p * Math.sqrt((2 * xsqr) / Pi);
    }

    let k = dof;

    while (k >= 2) {
      p = (p * xsqr) / k;
      k = k - 2;
    }

    let t = p;
    let a = dof;

    while (t > 0.0000000001 * p) {
      a = a + 2;
      t = (t * xsqr) / a;
      p = p + t;
    }

    return 1 - p
  }

  return Math.round(ChiSq(xsqr, dof) * 1000000) / 1000000
};

const CONFIDENCE = {};

/**
 * Returns the confidence interval for a population mean.
 *
 * Category: Statistical
 *
 * @param {*} alpha The significance level used to compute the confidence level. The confidence level equals 100*(1 - alpha)%, or in other words, an alpha of 0.05 indicates a 95 percent confidence level.
 * @param {*} standard_dev The population standard deviation for the data range and is assumed to be known.
 * @param {*} size The sample size.
 * @returns
 */
CONFIDENCE.NORM = (alpha, standard_dev, size) => {
  alpha = parseNumber(alpha);
  standard_dev = parseNumber(standard_dev);
  size = parseNumber(size);

  if (anyIsError(alpha, standard_dev, size)) {
    return value
  }

  return jstat.normalci(1, alpha, standard_dev, size)[1] - 1
};

/**
 * Returns the confidence interval for a population mean, using a Student's t distribution.
 *
 * Category: Statistical
 *
 * @param {*} alpha The significance level used to compute the confidence level. The confidence level equals 100*(1 - alpha)%, or in other words, an alpha of 0.05 indicates a 95 percent confidence level.
 * @param {*} standard_dev The population standard deviation for the data range and is assumed to be known.
 * @param {*} size The sample size.
 * @returns
 */
CONFIDENCE.T = (alpha, standard_dev, size) => {
  alpha = parseNumber(alpha);
  standard_dev = parseNumber(standard_dev);
  size = parseNumber(size);

  if (anyIsError(alpha, standard_dev, size)) {
    return value
  }

  return jstat.tci(1, alpha, standard_dev, size)[1] - 1
};

/**
 * Returns the correlation coefficient between two data sets.
 *
 * Category: Statistical
 *
 * @param {*} array1 A range of value values.
 * @param {*} array2 A second range of value values.
 * @returns
 */
function CORREL(array1, array2) {
  array1 = parseNumberArray(flatten(array1));
  array2 = parseNumberArray(flatten(array2));

  if (anyIsError(array1, array2)) {
    return value
  }

  return jstat.corrcoeff(array1, array2)
}

/**
 * Counts how many numbers are in the list of arguments.
 *
 * Category: Statistical
 *
 * @param {*} args Cell reference, or range within which you want to count numbers.count numbers.
 * @returns
 */
function COUNT() {
  const flatArguments = flatten(arguments);

  return numbers(flatArguments).length
}

/**
 * Counts how many values are in the list of arguments.
 *
 * Category: Statistical
 *
 * @param {*} args Arguments representing the values that you want to count.
 * @returns
 */
function COUNTA() {
  const flatArguments = flatten(arguments);

  return flatArguments.length - COUNTBLANK(flatArguments)
}

/**
 * Counts the number of blank values within a range.
 *
 * Category: Statistical
 *
 * @param {*} args The range from which you want to count the blank values.
 * @returns
 */
function COUNTBLANK() {
  const range = flatten(arguments);
  let blanks = 0;
  let element;

  for (let i = 0; i < range.length; i++) {
    element = range[i];

    if (element === undefined || element === null || element === '') {
      blanks++;
    }
  }

  return blanks
}

/**
 * Counts the number of values within a range that meet the given criteria.
 *
 * Category: Statistical
 *
 * @returns
 */
function COUNTIF(range, criteria) {
  range = flatten(range);

  const isWildcard = criteria === void 0 || criteria === '*';

  if (isWildcard) {
    return range.length
  }

  let matches = 0;
  const tokenizedCriteria = parse(criteria + '');

  for (let i = 0; i < range.length; i++) {
    const value = range[i];
    const tokens = [createToken(value, TOKEN_TYPE_LITERAL)].concat(tokenizedCriteria);

    if (compute(tokens)) {
      matches++;
    }
  }

  return matches
}

/**
 * Counts the number of values within a range that meet multiple criteria.
 *
 * Category: Statistical
 *
 * @param {*} args Range in which to evaluate the associated criteria.
 * @returns
 */
function COUNTIFS() {
  const args = argsToArray(arguments);
  const results = new Array(flatten(args[0]).length);

  for (let i = 0; i < results.length; i++) {
    results[i] = true;
  }

  for (let i = 0; i < args.length; i += 2) {
    const range = flatten(args[i]);
    const criteria = args[i + 1];
    const isWildcard = criteria === void 0 || criteria === '*';

    if (!isWildcard) {
      const tokenizedCriteria = parse(criteria + '');

      for (let j = 0; j < range.length; j++) {
        const value = range[j];
        const tokens = [createToken(value, TOKEN_TYPE_LITERAL)].concat(tokenizedCriteria);

        results[j] = results[j] && compute(tokens);
      }
    }
  }

  let result = 0;

  for (let i = 0; i < results.length; i++) {
    if (results[i]) {
      result++;
    }
  }

  return result
}

const COVARIANCE = {};

/**
 * Returns covariance, the average of the products of paired deviations.
 *
 * Category: Statistical
 *
 * @param {*} array1 The first value range of integers.
 * @param {*} array2 The second value range of integers.
 * @returns
 */
COVARIANCE.P = (array1, array2) => {
  array1 = parseNumberArray(flatten(array1));
  array2 = parseNumberArray(flatten(array2));

  if (anyIsError(array1, array2)) {
    return value
  }

  const mean1 = jstat.mean(array1);
  const mean2 = jstat.mean(array2);
  let result = 0;
  const n = array1.length;

  for (let i = 0; i < n; i++) {
    result += (array1[i] - mean1) * (array2[i] - mean2);
  }

  return result / n
};

/**
 * Returns the sample covariance, the average of the products deviations for each data point pair in two data sets.
 *
 * Category: Statistical
 *
 * @param {*} array1 The first value range of integers.
 * @param {*} array2 The second value range of integers.
 * @returns
 */
COVARIANCE.S = (array1, array2) => {
  array1 = parseNumberArray(flatten(array1));
  array2 = parseNumberArray(flatten(array2));

  if (anyIsError(array1, array2)) {
    return value
  }

  return jstat.covariance(array1, array2)
};

/**
 * Returns the sum of squares of deviations.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 arguments for which you want to calculate the sum of squared deviations. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
function DEVSQ() {
  const range = parseNumberArray(flatten(arguments));

  if (range instanceof Error) {
    return range
  }

  const mean = jstat.mean(range);
  let result = 0;

  for (let i = 0; i < range.length; i++) {
    result += Math.pow(range[i] - mean, 2);
  }

  return result
}

const EXPON = {};

/**
 * Returns the exponential distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The value of the function.
 * @param {*} lambda The parameter value.
 * @param {*} cumulative A logical value that indicates which form of the exponential function to provide. If cumulative is TRUE, EXPON.DIST returns the cumulative distribution function; if FALSE, it returns the probability density function.
 * @returns
 */
EXPON.DIST = (x, lambda, cumulative) => {
  x = parseNumber(x);
  lambda = parseNumber(lambda);

  if (anyIsError(x, lambda)) {
    return value
  }

  return cumulative ? jstat.exponential.cdf(x, lambda) : jstat.exponential.pdf(x, lambda)
};

const F = {};

/**
 * Returns the F probability distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The value at which to evaluate the function.
 * @param {*} deg_freedom1 The numerator degrees of freedom.
 * @param {*} deg_freedom2 The denominator degrees of freedom.
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, F.DIST returns the cumulative distribution function; if FALSE, it returns the probability density function.
 * @returns
 */
F.DIST = (x, deg_freedom1, deg_freedom2, cumulative) => {
  x = parseNumber(x);
  deg_freedom1 = parseNumber(deg_freedom1);
  deg_freedom2 = parseNumber(deg_freedom2);

  if (anyIsError(x, deg_freedom1, deg_freedom2)) {
    return value
  }

  return cumulative
    ? jstat.centralF.cdf(x, deg_freedom1, deg_freedom2)
    : jstat.centralF.pdf(x, deg_freedom1, deg_freedom2)
};

/**
 * Returns the F probability distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The value at which to evaluate the function.
 * @param {*} deg_freedom1 The numerator degrees of freedom.
 * @param {*} deg_freedom2 The denominator degrees of freedom.
 * @returns
 */
F.DIST.RT = function (x, deg_freedom1, deg_freedom2) {
  if (arguments.length !== 3) {
    return na
  }

  if (x < 0 || deg_freedom1 < 1 || deg_freedom2 < 1) {
    return num
  }

  if (typeof x !== 'number' || typeof deg_freedom1 !== 'number' || typeof deg_freedom2 !== 'number') {
    return value
  }

  return 1 - jstat.centralF.cdf(x, deg_freedom1, deg_freedom2)
};

/**
 * Returns the inverse of the F probability distribution.
 *
 * Category: Statistical
 *
 * @param {*} probability A probability associated with the F cumulative distribution.
 * @param {*} deg_freedom1 The numerator degrees of freedom.
 * @param {*} deg_freedom2 The denominator degrees of freedom.
 * @returns
 */
F.INV = (probability, deg_freedom1, deg_freedom2) => {
  probability = parseNumber(probability);
  deg_freedom1 = parseNumber(deg_freedom1);
  deg_freedom2 = parseNumber(deg_freedom2);

  if (anyIsError(probability, deg_freedom1, deg_freedom2)) {
    return value
  }

  if (probability <= 0.0 || probability > 1.0) {
    return num
  }

  return jstat.centralF.inv(probability, deg_freedom1, deg_freedom2)
};

/**
 * Returns the inverse of the F probability distribution.
 *
 * Category: Statistical
 *
 * @param {*} probability A probability associated with the F cumulative distribution.
 * @param {*} deg_freedom1 The numerator degrees of freedom.
 * @param {*} deg_freedom2 The denominator degrees of freedom.
 * @returns
 */
F.INV.RT = function (probability, deg_freedom1, deg_freedom2) {
  if (arguments.length !== 3) {
    return na
  }

  if (
    probability < 0 ||
    probability > 1 ||
    deg_freedom1 < 1 ||
    deg_freedom1 > Math.pow(10, 10) ||
    deg_freedom2 < 1 ||
    deg_freedom2 > Math.pow(10, 10)
  ) {
    return num
  }

  if (typeof probability !== 'number' || typeof deg_freedom1 !== 'number' || typeof deg_freedom2 !== 'number') {
    return value
  }

  return jstat.centralF.inv(1.0 - probability, deg_freedom1, deg_freedom2)
};

/**
 * Returns the result of an F-test.
 *
 * Category: Statistical
 *
 * @param {*} array1 The first array or range of data.
 * @param {*} array2 The second array or range of data.
 * @returns
 */
F.TEST = (array1, array2) => {
  if (!array1 || !array2) {
    return na
  }

  if (!(array1 instanceof Array) || !(array2 instanceof Array)) {
    return na
  }

  if (array1.length < 2 || array2.length < 2) {
    return div0
  }

  const sumOfSquares = (values, x1) => {
    let sum = 0;

    for (let i = 0; i < values.length; i++) {
      sum += Math.pow(values[i] - x1, 2);
    }

    return sum
  };

  const x1 = SUM(array1) / array1.length;
  const x2 = SUM(array2) / array2.length;
  const sum1 = sumOfSquares(array1, x1) / (array1.length - 1);
  const sum2 = sumOfSquares(array2, x2) / (array2.length - 1);

  return sum1 / sum2
};

/**
 * Returns the Fisher transformation.
 *
 * Category: Statistical
 *
 * @param {*} x A numeric value for which you want the transformation.
 * @returns
 */
function FISHER(x) {
  x = parseNumber(x);

  if (x instanceof Error) {
    return x
  }

  return Math.log((1 + x) / (1 - x)) / 2
}

/**
 * Returns the inverse of the Fisher transformation.
 *
 * Category: Statistical
 *
 * @param {*} y The value for which you want to perform the inverse of the transformation.
 * @returns
 */
function FISHERINV(y) {
  y = parseNumber(y);

  if (y instanceof Error) {
    return y
  }

  const e2y = Math.exp(2 * y);

  return (e2y - 1) / (e2y + 1)
}

/**
 * Returns a value along a linear trend.
 *
 * Category: Statistical
 *
 * @param {*} x The data point for which you want to predict a value.
 * @param {*} known_ys The dependent array or range of data.
 * @param {*} known_xs The independent array or range of data.
 * @returns
 */
function FORECAST(x, known_ys, known_xs) {
  x = parseNumber(x);
  known_ys = parseNumberArray(flatten(known_ys));
  known_xs = parseNumberArray(flatten(known_xs));

  if (anyIsError(x, known_ys, known_xs)) {
    return value
  }

  const xmean = jstat.mean(known_xs);
  const ymean = jstat.mean(known_ys);
  const n = known_xs.length;
  let num = 0;
  let den = 0;

  for (let i = 0; i < n; i++) {
    num += (known_xs[i] - xmean) * (known_ys[i] - ymean);
    den += Math.pow(known_xs[i] - xmean, 2);
  }

  const b = num / den;
  const a = ymean - b * xmean;

  return a + b * x
}

/**
 * Returns a frequency distribution as a vertical array.
 *
 * Category: Statistical
 *
 * @param {*} data_array An array of or reference to a set of values for which you want to count frequencies. If data_array contains no values, FREQUENCY returns an array of zeros.
 * @param {*} bins_array An array of or reference to intervals into which you want to group the values in data_array. If bins_array contains no values, FREQUENCY returns the number of elements in data_array.
 * @returns
 */
function FREQUENCY(data_array, bins_array) {
  data_array = parseNumberArray(flatten(data_array));
  bins_array = parseNumberArray(flatten(bins_array));

  if (anyIsError(data_array, bins_array)) {
    return value
  }

  const n = data_array.length;
  const b = bins_array.length;
  const r = [];

  for (let i = 0; i <= b; i++) {
    r[i] = 0;

    for (let j = 0; j < n; j++) {
      if (i === 0) {
        if (data_array[j] <= bins_array[0]) {
          r[0] += 1;
        }
      } else if (i < b) {
        if (data_array[j] > bins_array[i - 1] && data_array[j] <= bins_array[i]) {
          r[i] += 1;
        }
      } else if (i === b) {
        if (data_array[j] > bins_array[b - 1]) {
          r[b] += 1;
        }
      }
    }
  }

  return r
}

/**
 * Returns the Gamma function value.
 *
 * Category: Statistical
 *
 * @param {*} number Returns a number.
 * @returns
 */
function GAMMA(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  if (number === 0) {
    return num
  }

  if (parseInt(number, 10) === number && number < 0) {
    return num
  }

  return jstat.gammafn(number)
}

/**
 * Returns the gamma distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The value at which you want to evaluate the distribution.
 * @param {*} alpha A parameter to the distribution.
 * @param {*} beta A parameter to the distribution. If beta = 1, GAMMA.DIST returns the standard gamma distribution.
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, GAMMA.DIST returns the cumulative distribution function; if FALSE, it returns the probability density function.
 * @returns
 */
GAMMA.DIST = function (value$1, alpha, beta, cumulative) {
  if (arguments.length !== 4) {
    return na
  }

  if (value$1 < 0 || alpha <= 0 || beta <= 0) {
    return value
  }

  if (typeof value$1 !== 'number' || typeof alpha !== 'number' || typeof beta !== 'number') {
    return value
  }

  return cumulative ? jstat.gamma.cdf(value$1, alpha, beta, true) : jstat.gamma.pdf(value$1, alpha, beta, false)
};

/**
 * Returns the inverse of the gamma cumulative distribution.
 *
 * Category: Statistical
 *
 * @param {*} probability The probability associated with the gamma distribution.
 * @param {*} alpha A parameter to the distribution.
 * @param {*} beta A parameter to the distribution. If beta = 1, GAMMA.INV returns the standard gamma distribution.
 * @returns
 */
GAMMA.INV = function (probability, alpha, beta) {
  if (arguments.length !== 3) {
    return na
  }

  if (probability < 0 || probability > 1 || alpha <= 0 || beta <= 0) {
    return num
  }

  if (typeof probability !== 'number' || typeof alpha !== 'number' || typeof beta !== 'number') {
    return value
  }

  return jstat.gamma.inv(probability, alpha, beta)
};

/**
 * Returns the natural logarithm of the gamma function, Γ(x).
 *
 * Category: Statistical
 *
 * @param {*} x The value for which you want to calculate GAMMALN.
 * @returns
 */
function GAMMALN(x) {
  x = parseNumber(x);

  if (x instanceof Error) {
    return x
  }

  return jstat.gammaln(x)
}

/**
 * Returns the natural logarithm of the gamma function, Γ(x).
 *
 * Category: Statistical
 *
 * @param {*} x The value for which you want to calculate GAMMALN.PRECISE.
 * @returns
 */
GAMMALN.PRECISE = function (x) {
  if (arguments.length !== 1) {
    return na
  }

  if (x <= 0) {
    return num
  }

  if (typeof x !== 'number') {
    return value
  }

  return jstat.gammaln(x)
};

/**
 * Returns 0.5 less than the standard normal cumulative distribution.
 *
 * Category: Statistical
 *
 * @param {*} z Returns a number.
 * @returns
 */
function GAUSS(z) {
  z = parseNumber(z);

  if (z instanceof Error) {
    return z
  }

  return jstat.normal.cdf(z, 0, 1) - 0.5
}

/**
 * Returns the geometric mean.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 arguments for which you want to calculate the mean. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
function GEOMEAN() {
  const args = parseNumberArray(flatten(arguments));

  if (args instanceof Error) {
    return args
  }

  return jstat.geomean(args)
}

/**
 * Returns values along an exponential trend.
 *
 * Category: Statistical
 *
 * @param {*} known_y The set of y-values you already know in the relationship y = b*m^x.
 - If the array known_y's is in a single column, then each column of known_x's is interpreted as a separate variable.
 - If the array known_y's is in a single row, then each row of known_x's is interpreted as a separate variable.
 - If any of the numbers in known_y's is 0 or negative, GROWTH returns the #NUM! error value.
 * @param {*} known_x Optional. An optional set of x-values that you may already know in the relationship y = b*m^x.
 - The array known_x's can include one or more sets of variables. If only one variable is used, known_y's and known_x's can be ranges of any shape, as long as they have equal dimensions. If more than one variable is used, known_y's must be a vector (that is, a range with a height of one row or a width of one column).
 - If known_x's is omitted, it is assumed to be the array {1,2,3,...} that is the same size as known_y's.
 * @param {*} new_x Optional. Are new x-values for which you want GROWTH to return corresponding y-values.
 - new_x's must include a column (or row) for each independent variable, just as known_x's does. So, if known_y's is in a single column, known_x's and new_x's must have the same number of columns. If known_y's is in a single row, known_x's and new_x's must have the same number of rows.
 - If new_x's is omitted, it is assumed to be the same as known_x's.
 - If both known_x's and new_x's are omitted, they are assumed to be the array {1,2,3,...} that is the same size as known_y's.
 * @param {*} use_const Optional. A logical value specifying whether to force the constant b to equal 1. If const is TRUE or omitted, b is calculated normally. If const is FALSE, b is set equal to 1 and the m-values are adjusted so that y = m^x.
 - If const is TRUE or omitted, b is calculated normally.
 - If const is FALSE, b is set equal to 1 and the m-values are adjusted so that y = m^x.
 * @returns
 */
function GROWTH(known_y, known_x, new_x, use_const) {
  // Credits: Ilmari Karonen (http://stackoverflow.com/questions/14161990/how-to-implement-growth-function-in-javascript)
  known_y = parseNumberArray(known_y);

  if (known_y instanceof Error) {
    return known_y
  }

  // Default values for optional parameters:
  let i;

  if (known_x === undefined) {
    known_x = [];

    for (i = 1; i <= known_y.length; i++) {
      known_x.push(i);
    }
  }

  if (new_x === undefined) {
    new_x = [];

    for (i = 1; i <= known_y.length; i++) {
      new_x.push(i);
    }
  }

  known_x = parseNumberArray(known_x);
  new_x = parseNumberArray(new_x);

  if (anyIsError(known_x, new_x)) {
    return value
  }

  if (use_const === undefined) {
    use_const = true;
  }

  // Calculate sums over the data:
  const n = known_y.length;
  let avg_x = 0;
  let avg_y = 0;
  let avg_xy = 0;
  let avg_xx = 0;

  for (i = 0; i < n; i++) {
    const x = known_x[i];
    const y = Math.log(known_y[i]);
    avg_x += x;
    avg_y += y;
    avg_xy += x * y;
    avg_xx += x * x;
  }

  avg_x /= n;
  avg_y /= n;
  avg_xy /= n;
  avg_xx /= n;

  // Compute linear regression coefficients:
  let beta;
  let alpha;

  if (use_const) {
    beta = (avg_xy - avg_x * avg_y) / (avg_xx - avg_x * avg_x);
    alpha = avg_y - beta * avg_x;
  } else {
    beta = avg_xy / avg_xx;
    alpha = 0;
  }

  // Compute and return result array:
  const new_y = [];

  for (i = 0; i < new_x.length; i++) {
    new_y.push(Math.exp(alpha + beta * new_x[i]));
  }

  return new_y
}

/**
 * Returns the harmonic mean.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 arguments for which you want to calculate the mean. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
function HARMEAN() {
  const range = parseNumberArray(flatten(arguments));

  if (range instanceof Error) {
    return range
  }

  const n = range.length;
  let den = 0;

  for (let i = 0; i < n; i++) {
    den += 1 / range[i];
  }

  return n / den
}

const HYPGEOM = {};

/**
 * Returns the hypergeometric distribution.
 *
 * Category: Statistical
 *
 * @param {*} sample_s The number of successes in the sample.
 * @param {*} number_sample The size of the sample.
 * @param {*} population_s The number of successes in the population.
 * @param {*} number_pop The population size.
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, then HYPGEOM.DIST returns the cumulative distribution function; if FALSE, it returns the probability mass function.
 * @returns
 */
HYPGEOM.DIST = (sample_s, number_sample, population_s, number_pop, cumulative) => {
  sample_s = parseNumber(sample_s);
  number_sample = parseNumber(number_sample);
  population_s = parseNumber(population_s);
  number_pop = parseNumber(number_pop);

  if (anyIsError(sample_s, number_sample, population_s, number_pop)) {
    return value
  }

  function pdf(x, n, M, N) {
    return (COMBIN(M, x) * COMBIN(N - M, n - x)) / COMBIN(N, n)
  }

  function cdf(x, n, M, N) {
    let result = 0;

    for (let i = 0; i <= x; i++) {
      result += pdf(i, n, M, N);
    }

    return result
  }

  return cumulative
    ? cdf(sample_s, number_sample, population_s, number_pop)
    : pdf(sample_s, number_sample, population_s, number_pop)
};

/**
 * Returns the intercept of the linear regression line.
 *
 * Category: Statistical
 *
 * @param {*} known_y The dependent set of observations or data.
 * @param {*} known_x The independent set of observations or data.
 * @returns
 */
function INTERCEPT(known_y, known_x) {
  known_y = parseNumberArray(known_y);
  known_x = parseNumberArray(known_x);

  if (anyIsError(known_y, known_x)) {
    return value
  }

  if (known_y.length !== known_x.length) {
    return na
  }

  return FORECAST(0, known_y, known_x)
}

/**
 * Returns the kurtosis of a data set.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 arguments for which you want to calculate kurtosis. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
function KURT() {
  const range = parseNumberArray(flatten(arguments));

  if (range instanceof Error) {
    return range
  }

  const mean = jstat.mean(range);
  const n = range.length;
  let sigma = 0;

  for (let i = 0; i < n; i++) {
    sigma += Math.pow(range[i] - mean, 4);
  }

  sigma = sigma / Math.pow(jstat.stdev(range, true), 4);

  return ((n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))) * sigma - (3 * (n - 1) * (n - 1)) / ((n - 2) * (n - 3))
}

/**
 * Returns the k-th largest value in a data set.
 *
 * Category: Statistical
 *
 * @param {*} array The array or range of data for which you want to determine the k-th largest value.
 * @param {*} k The position (from the largest) in the array or value range of data to return.
 * @returns
 */
function LARGE(array, k) {
  array = parseNumberArray(flatten(array));
  k = parseNumber(k);

  if (anyIsError(array, k)) {
    return array
  }

  if (k < 0 || array.length < k) {
    return value
  }

  return array.sort((a, b) => b - a)[k - 1]
}

/**
 * Returns the parameters of a linear trend.
 *
 * Category: Statistical
 *
 * @param {*} known_y The set of y-values that you already know in the relationship y = mx + b.
 - If the range of known_y's is in a single column, each column of known_x's is interpreted as a separate variable.
 - If the range of known_y's is contained in a single row, each row of known_x's is interpreted as a separate variable.
 * @param {*} known_x Optional. A set of x-values that you may already know in the relationship y = mx + b.
 - The range of known_x's can include one or more sets of variables. If only one variable is used, known_y's and known_x's can be ranges of any shape, as long as they have equal dimensions. If more than one variable is used, known_y's must be a vector (that is, a range with a height of one row or a width of one column).
 - If known_x's is omitted, it is assumed to be the array {1,2,3,...} that is the same size as known_y's.
 * @returns
 */
function LINEST(known_y, known_x) {
  known_y = parseNumberArray(flatten(known_y));
  known_x = parseNumberArray(flatten(known_x));

  if (anyIsError(known_y, known_x)) {
    return value
  }

  const ymean = jstat.mean(known_y);
  const xmean = jstat.mean(known_x);
  const n = known_x.length;
  let num = 0;
  let den = 0;

  for (let i = 0; i < n; i++) {
    num += (known_x[i] - xmean) * (known_y[i] - ymean);
    den += Math.pow(known_x[i] - xmean, 2);
  }

  const m = num / den;
  const b = ymean - m * xmean;

  return [m, b]
}

// According to Microsoft:
// http://office.microsoft.com/en-us/starter-help/logest-function-HP010342665.aspx
// LOGEST returns are based on the following linear model:
// ln y = x1 ln m1 + ... + xn ln mn + ln b
/**
 * Returns the parameters of an exponential trend.
 *
 * Category: Statistical
 *
 * @param {*} known_y The set of y-values you already know in the relationship y = b*m^x.
 - If the array known_y's is in a single column, then each column of known_x's is interpreted as a separate variable.
 - If the array known_y's is in a single row, then each row of known_x's is interpreted as a separate variable.
 * @param {*} known_x Optional. An optional set of x-values that you may already know in the relationship y = b*m^x.
 - The array known_x's can include one or more sets of variables. If only one variable is used, known_y's and known_x's can be ranges of any shape, as long as they have equal dimensions. If more than one variable is used, known_y's must be a range of values with a height of one row or a width of one column (which is also known as a vector).
 - If known_x's is omitted, it is assumed to be the array {1,2,3,...} that is the same size as known_y's.
 * @returns
 */
function LOGEST(known_y, known_x) {
  known_y = parseNumberArray(flatten(known_y));
  known_x = parseNumberArray(flatten(known_x));

  if (anyIsError(known_y, known_x)) {
    return value
  }

  if (known_y.length !== known_x.length) {
    return value
  }

  for (let i = 0; i < known_y.length; i++) {
    known_y[i] = Math.log(known_y[i]);
  }

  const result = LINEST(known_y, known_x);
  result[0] = Math.round(Math.exp(result[0]) * 1000000) / 1000000;
  result[1] = Math.round(Math.exp(result[1]) * 1000000) / 1000000;

  return result
}

const LOGNORM = {};

/**
 * Returns the cumulative lognormal distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The value at which to evaluate the function.
 * @param {*} mean The mean of ln(x).
 * @param {*} standard_dev The standard deviation of ln(x).
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, LOGNORM.DIST returns the cumulative distribution function; if FALSE, it returns the probability density function.
 * @returns
 */
LOGNORM.DIST = (x, mean, standard_dev, cumulative) => {
  x = parseNumber(x);
  mean = parseNumber(mean);
  standard_dev = parseNumber(standard_dev);

  if (anyIsError(x, mean, standard_dev)) {
    return value
  }

  return cumulative ? jstat.lognormal.cdf(x, mean, standard_dev) : jstat.lognormal.pdf(x, mean, standard_dev)
};

/**
 * Returns the inverse of the lognormal cumulative distribution.
 *
 * Category: Statistical
 *
 * @param {*} probability A probability associated with the lognormal distribution.
 * @param {*} mean The mean of ln(x).
 * @param {*} standard_dev The standard deviation of ln(x).
 * @returns
 */
LOGNORM.INV = (probability, mean, standard_dev) => {
  probability = parseNumber(probability);
  mean = parseNumber(mean);
  standard_dev = parseNumber(standard_dev);

  if (anyIsError(probability, mean, standard_dev)) {
    return value
  }

  return jstat.lognormal.inv(probability, mean, standard_dev)
};

/**
 * Returns the maximum value in a list of arguments.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 numbers for which you want to find the maximum value.
 * @returns
 */
function MAX() {
  const flatArguments = flatten(arguments);
  const someError = anyError.apply(undefined, flatArguments);

  if (someError) {
    return someError
  }

  const range = numbers(flatArguments);

  return range.length === 0 ? 0 : Math.max.apply(Math, range)
}

/**
 * Returns the maximum value in a list of arguments, including numbers, text, and logical values.
 *
 * Category: Statistical
 *
 * @param {*} args value1, value2,... Number arguments 2 to 255 for which you want to find the largest value.
 * @returns
 */
function MAXA() {
  const flatArguments = flatten(arguments);
  const someError = anyError.apply(undefined, flatArguments);

  if (someError) {
    return someError
  }

  let range = arrayValuesToNumbers(flatArguments);
  range = range.map((value) => (value === undefined || value === null ? 0 : value));

  return range.length === 0 ? 0 : Math.max.apply(Math, range)
}

/**
 * Returns the median of the given numbers.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 numbers for which you want the median.
 * @returns
 */
function MEDIAN() {
  const flatArguments = flatten(arguments);
  const someError = anyError.apply(undefined, flatArguments);

  if (someError) {
    return someError
  }

  const range = arrayValuesToNumbers(flatArguments);
  let result = jstat.median(range);

  if (isNaN(result)) {
    result = num;
  }

  return result
}

/**
 * Returns the minimum value in a list of arguments.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number1 is optional, subsequent numbers are optional. 1 to 255 numbers for which you want to find the minimum value.
 * @returns
 */
function MIN() {
  const flatArguments = flatten(arguments);
  const someError = anyError.apply(undefined, flatArguments);

  if (someError) {
    return someError
  }

  const range = numbers(flatArguments);

  return range.length === 0 ? 0 : Math.min.apply(Math, range)
}

/**
 * Returns the smallest value in a list of arguments, including numbers, text, and logical values.
 *
 * Category: Statistical
 *
 * @param {*} args value1, value2, ... Value1 is required, subsequent values are optional. 1 to 255 values for which you want to find the smallest value.
 * @returns
 */
function MINA() {
  const flatArguments = flatten(arguments);
  const someError = anyError.apply(undefined, flatArguments);

  if (someError) {
    return someError
  }

  let range = arrayValuesToNumbers(flatArguments);
  range = range.map((value) => (value === undefined || value === null ? 0 : value));

  return range.length === 0 ? 0 : Math.min.apply(Math, range)
}

const MODE = {};

/**
 * Returns a vertical array of the most frequently occurring, or repetitive values in an array or range of data.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number arguments 2 to 254 for which you want to calculate the mode. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
MODE.MULT = function () {
  // Credits: Roönaän
  const range = parseNumberArray(flatten(arguments));

  if (range instanceof Error) {
    return range
  }

  const n = range.length;
  const count = {};
  let maxItems = [];
  let max = 0;
  let currentItem;

  for (let i = 0; i < n; i++) {
    currentItem = range[i];
    count[currentItem] = count[currentItem] ? count[currentItem] + 1 : 1;

    if (count[currentItem] > max) {
      max = count[currentItem];
      maxItems = [];
    }

    if (count[currentItem] === max) {
      maxItems[maxItems.length] = currentItem;
    }
  }

  return maxItems
};

/**
 * Returns the most common value in a data set.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Arguments 2 to 254 for which you want to calculate the mode. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
MODE.SNGL = function () {
  const range = parseNumberArray(flatten(arguments));

  if (range instanceof Error) {
    return range
  }

  return MODE.MULT(range).sort((a, b) => a - b)[0]
};

const NEGBINOM = {};

/**
 * Returns the negative binomial distribution.
 *
 * Category: Statistical
 *
 * @param {*} number_f The number of failures.
 * @param {*} number_s The threshold number of successes.
 * @param {*} probability_s The probability of a success.
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, NEGBINOM.DIST returns the cumulative distribution function; if FALSE, it returns the probability density function.
 * @returns
 */
NEGBINOM.DIST = (number_f, number_s, probability_s, cumulative) => {
  number_f = parseNumber(number_f);
  number_s = parseNumber(number_s);
  probability_s = parseNumber(probability_s);

  if (anyIsError(number_f, number_s, probability_s)) {
    return value
  }

  return cumulative
    ? jstat.negbin.cdf(number_f, number_s, probability_s)
    : jstat.negbin.pdf(number_f, number_s, probability_s)
};

const NORM = {};

/**
 * Returns the normal cumulative distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The value for which you want the distribution.
 * @param {*} mean The arithmetic mean of the distribution.
 * @param {*} standard_dev The standard deviation of the distribution.
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, NORM.DIST returns the cumulative distribution function; if FALSE, it returns the probability density function.
 * @returns
 */
NORM.DIST = (x, mean, standard_dev, cumulative) => {
  x = parseNumber(x);
  mean = parseNumber(mean);
  standard_dev = parseNumber(standard_dev);

  if (anyIsError(x, mean, standard_dev)) {
    return value
  }

  if (standard_dev <= 0) {
    return num
  }

  // Return normal distribution computed by jStat [http://jstat.org]
  return cumulative ? jstat.normal.cdf(x, mean, standard_dev) : jstat.normal.pdf(x, mean, standard_dev)
};

/**
 * Returns the inverse of the normal cumulative distribution.
 *
 * Category: Statistical
 *
 * @param {*} probability A probability corresponding to the normal distribution.
 * @param {*} mean The arithmetic mean of the distribution.
 * @param {*} standard_dev The standard deviation of the distribution.
 * @returns
 */
NORM.INV = (probability, mean, standard_dev) => {
  probability = parseNumber(probability);
  mean = parseNumber(mean);
  standard_dev = parseNumber(standard_dev);

  if (anyIsError(probability, mean, standard_dev)) {
    return value
  }

  return jstat.normal.inv(probability, mean, standard_dev)
};

NORM.S = {};

/**
 * Returns the standard normal cumulative distribution.
 *
 * Category: Statistical
 *
 * @param {*} z The value for which you want the distribution.
 * @param {*} cumulative Cumulative is a logical value that determines the form of the function. If cumulative is TRUE, NORMS.DIST returns the cumulative distribution function; if FALSE, it returns the probability mass function.
 * @returns
 */
NORM.S.DIST = (z, cumulative) => {
  z = parseNumber(z);

  if (z instanceof Error) {
    return value
  }

  return cumulative ? jstat.normal.cdf(z, 0, 1) : jstat.normal.pdf(z, 0, 1)
};

/**
 * Returns the inverse of the standard normal cumulative distribution.
 *
 * Category: Statistical
 *
 * @param {*} probability A probability corresponding to the normal distribution.
 * @returns
 */
NORM.S.INV = (probability) => {
  probability = parseNumber(probability);

  if (probability instanceof Error) {
    return value
  }

  return jstat.normal.inv(probability, 0, 1)
};

/**
 * Returns the Pearson product moment correlation coefficient.
 *
 * Category: Statistical
 *
 * @param {*} array1 A set of independent values.
 * @param {*} array2 A set of dependent values.
 * @returns
 */
function PEARSON(array1, array2) {
  array2 = parseNumberArray(flatten(array2));
  array1 = parseNumberArray(flatten(array1));

  if (anyIsError(array2, array1)) {
    return value
  }

  const xmean = jstat.mean(array1);
  const ymean = jstat.mean(array2);
  const n = array1.length;
  let num = 0;
  let den1 = 0;
  let den2 = 0;

  for (let i = 0; i < n; i++) {
    num += (array1[i] - xmean) * (array2[i] - ymean);
    den1 += Math.pow(array1[i] - xmean, 2);
    den2 += Math.pow(array2[i] - ymean, 2);
  }

  return num / Math.sqrt(den1 * den2)
}

const PERCENTILE = {};

/**
 * Returns the k-th percentile of values in a range, where k is in the range 0..1, exclusive.
 *
 * Category: Statistical
 *
 * @returns
 */
PERCENTILE.EXC = (array, k) => {
  array = parseNumberArray(flatten(array));
  k = parseNumber(k);

  if (anyIsError(array, k)) {
    return value
  }

  array = array.sort((a, b) => a - b);
  const n = array.length;

  if (k < 1 / (n + 1) || k > 1 - 1 / (n + 1)) {
    return num
  }

  const l = k * (n + 1) - 1;
  const fl = Math.floor(l);

  return cleanFloat(l === fl ? array[l] : array[fl] + (l - fl) * (array[fl + 1] - array[fl]))
};

/**
 * Returns the k-th percentile of values in a range.
 *
 * Category: Statistical
 *
 * @param {*} array The array or range of data that defines relative standing.
 * @param {*} k The percentile value in the range 0..1, inclusive.
 * @returns
 */
PERCENTILE.INC = (array, k) => {
  array = parseNumberArray(flatten(array));
  k = parseNumber(k);

  if (anyIsError(array, k)) {
    return value
  }

  array = array.sort((a, b) => a - b);
  const n = array.length;
  const l = k * (n - 1);
  const fl = Math.floor(l);

  return cleanFloat(l === fl ? array[l] : array[fl] + (l - fl) * (array[fl + 1] - array[fl]))
};

const PERCENTRANK = {};

/**
 * Returns the rank of a value in a data set as a percentage (0..1, exclusive) of the data set.
 *
 * Category: Statistical
 *
 * @param {*} array The array or range of data with numeric values that defines relative standing
 * @param {*} x The value for which you want to know the rank.
 * @param {*} significance Optional. A value that identifies the number of significant digits for the returned percentage value. If omitted, PERCENTRANK.EXC uses three digits (0.xxx).
 * @returns
 */
PERCENTRANK.EXC = (array, x, significance) => {
  significance = significance === undefined ? 3 : significance;
  array = parseNumberArray(flatten(array));
  x = parseNumber(x);
  significance = parseNumber(significance);

  if (anyIsError(array, x, significance)) {
    return value
  }

  array = array.sort((a, b) => a - b);
  const uniques = UNIQUE.apply(null, array);
  const n = array.length;
  const m = uniques.length;
  const power = Math.pow(10, significance);
  let result = 0;
  let match = false;
  let i = 0;

  while (!match && i < m) {
    if (x === uniques[i]) {
      result = (array.indexOf(uniques[i]) + 1) / (n + 1);
      match = true;
    } else if (x >= uniques[i] && (x < uniques[i + 1] || i === m - 1)) {
      result = (array.indexOf(uniques[i]) + 1 + (x - uniques[i]) / (uniques[i + 1] - uniques[i])) / (n + 1);
      match = true;
    }

    i++;
  }

  return Math.floor(result * power) / power
};

/**
 * Returns the percentage rank of a value in a data set.
 *
 * Category: Statistical
 *
 * @param {*} array The array or range of data with numeric values that defines relative standing.
 * @param {*} x The value for which you want to know the rank.
 * @param {*} significance Optional. A value that identifies the number of significant digits for the returned percentage value. If omitted, PERCENTRANK.INC uses three digits (0.xxx).
 * @returns
 */
PERCENTRANK.INC = (array, x, significance) => {
  significance = significance === undefined ? 3 : significance;
  array = parseNumberArray(flatten(array));
  x = parseNumber(x);
  significance = parseNumber(significance);

  if (anyIsError(array, x, significance)) {
    return value
  }

  array = array.sort((a, b) => a - b);
  const uniques = UNIQUE.apply(null, array);
  const n = array.length;
  const m = uniques.length;
  const power = Math.pow(10, significance);
  let result = 0;
  let match = false;
  let i = 0;

  while (!match && i < m) {
    if (x === uniques[i]) {
      result = array.indexOf(uniques[i]) / (n - 1);
      match = true;
    } else if (x >= uniques[i] && (x < uniques[i + 1] || i === m - 1)) {
      result = (array.indexOf(uniques[i]) + (x - uniques[i]) / (uniques[i + 1] - uniques[i])) / (n - 1);
      match = true;
    }

    i++;
  }

  return Math.floor(result * power) / power
};

/**
 * Returns the number of permutations for a given number of objects.
 *
 * Category: Statistical
 *
 * @param {*} number An integer that describes the number of objects.
 * @param {*} number_chosen An integer that describes the number of objects in each permutation.
 * @returns
 */
function PERMUT(number, number_chosen) {
  number = parseNumber(number);
  number_chosen = parseNumber(number_chosen);

  if (anyIsError(number, number_chosen)) {
    return value
  }

  return FACT(number) / FACT(number - number_chosen)
}

/**
 * Returns the number of permutations for a given number of objects (with repetitions) that can be selected from the total objects.
 *
 * Category: Statistical
 *
 * @param {*} number An integer that describes the total number of objects.
 * @param {*} number_chosen An integer that describes the number of objects in each permutation.
 * @returns
 */
function PERMUTATIONA(number, number_chosen) {
  number = parseNumber(number);
  number_chosen = parseNumber(number_chosen);

  if (anyIsError(number, number_chosen)) {
    return value
  }

  return Math.pow(number, number_chosen)
}

/**
 * Returns the value of the density function for a standard normal distribution.
 *
 * Category: Statistical
 *
 * @param {*} x X is the number for which you want the density of the standard normal distribution.
 * @returns
 */
function PHI(x) {
  x = parseNumber(x);

  if (x instanceof Error) {
    return value
  }

  return Math.exp(-0.5 * x * x) / SQRT2PI
}

const POISSON = {};

/**
 * Returns the Poisson distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The number of events.
 * @param {*} mean The expected numeric value.
 * @param {*} cumulative A logical value that determines the form of the probability distribution returned. If cumulative is TRUE, POISSON.DIST returns the cumulative Poisson probability that the number of random events occurring will be between zero and x inclusive; if FALSE, it returns the Poisson probability mass function that the number of events occurring will be exactly x.
 * @returns
 */
POISSON.DIST = (x, mean, cumulative) => {
  x = parseNumber(x);
  mean = parseNumber(mean);

  if (anyIsError(x, mean)) {
    return value
  }

  return cumulative ? jstat.poisson.cdf(x, mean) : jstat.poisson.pdf(x, mean)
};

/**
 * Returns the probability that values in a range are between two limits.
 *
 * Category: Statistical
 *
 * @param {*} x_range The range of numeric values of x with which there are associated probabilities.
 * @param {*} prob_range A set of probabilities associated with values in x_range.
 * @param {*} lower_limit Optional. The lower bound on the value for which you want a probability.
 * @param {*} upper_limit Optional. The optional upper bound on the value for which you want a probability.
 * @returns
 */
function PROB(x_range, prob_range, lower_limit, upper_limit) {
  if (lower_limit === undefined) {
    return 0
  }

  upper_limit = upper_limit === undefined ? lower_limit : upper_limit;

  x_range = parseNumberArray(flatten(x_range));
  prob_range = parseNumberArray(flatten(prob_range));
  lower_limit = parseNumber(lower_limit);
  upper_limit = parseNumber(upper_limit);

  if (anyIsError(x_range, prob_range, lower_limit, upper_limit)) {
    return value
  }

  if (lower_limit === upper_limit) {
    return x_range.indexOf(lower_limit) >= 0 ? prob_range[x_range.indexOf(lower_limit)] : 0
  }

  const sorted = x_range.sort((a, b) => a - b);
  const n = sorted.length;
  let result = 0;

  for (let i = 0; i < n; i++) {
    if (sorted[i] >= lower_limit && sorted[i] <= upper_limit) {
      result += prob_range[x_range.indexOf(sorted[i])];
    }
  }

  return result
}

const QUARTILE = {};

/**
 * Returns the quartile of the data set, based on percentile values from 0..1, exclusive.
 *
 * Category: Statistical
 *
 * @param {*} array The array or value range of numeric values for which you want the quartile value.
 * @param {*} quart Indicates which value to return.
 * @returns
 */
QUARTILE.EXC = (range, quart) => {
  range = parseNumberArray(flatten(range));
  quart = parseNumber(quart);

  if (anyIsError(range, quart)) {
    return value
  }

  switch (quart) {
    case 1:
      return PERCENTILE.EXC(range, 0.25)
    case 2:
      return PERCENTILE.EXC(range, 0.5)
    case 3:
      return PERCENTILE.EXC(range, 0.75)
    default:
      return num
  }
};

/**
 * Returns the quartile of a data set.
 *
 * Category: Statistical
 *
 * @param {*} array The array or value range of numeric values for which you want the quartile value.
 * @param {*} quart Indicates which value to return.
 * @returns
 */
QUARTILE.INC = (range, quart) => {
  range = parseNumberArray(flatten(range));
  quart = parseNumber(quart);

  if (anyIsError(range, quart)) {
    return value
  }

  switch (quart) {
    case 1:
      return PERCENTILE.INC(range, 0.25)
    case 2:
      return PERCENTILE.INC(range, 0.5)
    case 3:
      return PERCENTILE.INC(range, 0.75)
    default:
      return num
  }
};

const RANK = {};

/**
 * Returns the rank of a number in a list of numbers.
 *
 * Category: Statistical
 *
 * @param {*} number The number whose rank you want to find.
 * @param {*} ref An array of, or a reference to, a list of numbers. Nonnumeric values in Ref are ignored.
 * @param {*} order Optional. A number specifying how to rank number.
 * @returns
 */
RANK.AVG = (number, ref, order) => {
  number = parseNumber(number);
  ref = parseNumberArray(flatten(ref));

  if (anyIsError(number, ref)) {
    return value
  }

  ref = flatten(ref);
  order = order || false;
  const sort = order ? (a, b) => a - b : (a, b) => b - a;
  ref = ref.sort(sort);

  const length = ref.length;
  let count = 0;

  for (let i = 0; i < length; i++) {
    if (ref[i] === number) {
      count++;
    }
  }

  return count > 1 ? (2 * ref.indexOf(number) + count + 1) / 2 : ref.indexOf(number) + 1
};

/**
 * Returns the rank of a number in a list of numbers.
 *
 * Category: Statistical
 *
 * @param {*} number The number whose rank you want to find.
 * @param {*} ref An array of, or a reference to, a list of numbers. Non-numeric values in Ref are ignored.
 * @param {*} order Optional. A number specifying how to rank number.
 * @returns
 */
RANK.EQ = (number, ref, order) => {
  number = parseNumber(number);
  ref = parseNumberArray(flatten(ref));

  if (anyIsError(number, ref)) {
    return value
  }

  order = order || false;
  const sort = order ? (a, b) => a - b : (a, b) => b - a;
  ref = ref.sort(sort);

  return ref.indexOf(number) + 1
};

/**
 * Returns the row number of a reference.
 *
 * Category: Lookup and reference
 *
 * @param {*} reference the value or range of values for which you want the row number.
 * @param {*} index
 * @returns
 */
function ROW(reference, index) {
  if (arguments.length !== 2) {
    return na
  }

  if (index < 0) {
    return num
  }

  if (!(reference instanceof Array) || typeof index !== 'number') {
    return value
  }

  if (reference.length === 0) {
    return undefined
  }

  return jstat.row(reference, index)
}

/**
 * Returns the square of the Pearson product moment correlation coefficient.
 *
 * Category: Statistical
 *
 * @param {*} known_y An array or range of data points.
 * @param {*} known_x An array or range of data points.
 * @returns
 */
function RSQ(known_y, known_x) {
  // no need to flatten here, PEARSON will take care of that
  known_y = parseNumberArray(flatten(known_y));
  known_x = parseNumberArray(flatten(known_x));

  if (anyIsError(known_y, known_x)) {
    return value
  }

  return Math.pow(PEARSON(known_y, known_x), 2)
}

/**
 * Returns the skewness of a distribution.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 arguments for which you want to calculate skewness. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
function SKEW() {
  const range = parseNumberArray(flatten(arguments));

  if (range instanceof Error) {
    return range
  }

  const mean = jstat.mean(range);
  const n = range.length;
  let sigma = 0;

  for (let i = 0; i < n; i++) {
    sigma += Math.pow(range[i] - mean, 3);
  }

  return (n * sigma) / ((n - 1) * (n - 2) * Math.pow(jstat.stdev(range, true), 3))
}

/**
 * Returns the skewness of a distribution based on a population.
 *
 * Category: Statistical
 *
 * @returns
 */
SKEW.P = function () {
  const range = parseNumberArray(flatten(arguments));

  if (range instanceof Error) {
    return range
  }

  const mean = jstat.mean(range);
  const n = range.length;
  let m2 = 0;
  let m3 = 0;

  for (let i = 0; i < n; i++) {
    m3 += Math.pow(range[i] - mean, 3);
    m2 += Math.pow(range[i] - mean, 2);
  }

  m3 = m3 / n;
  m2 = m2 / n;

  return m3 / Math.pow(m2, 3 / 2)
};

/**
 * Returns the slope of the linear regression line.
 *
 * Category: Statistical
 *
 * @param {*} known_y An array or value range of numeric dependent data points.
 * @param {*} known_x The set of independent data points.
 * @returns
 */
function SLOPE(known_y, known_x) {
  known_y = parseNumberArray(flatten(known_y));
  known_x = parseNumberArray(flatten(known_x));

  if (anyIsError(known_y, known_x)) {
    return value
  }

  const xmean = jstat.mean(known_x);
  const ymean = jstat.mean(known_y);
  const n = known_x.length;
  let num = 0;
  let den = 0;

  for (let i = 0; i < n; i++) {
    num += (known_x[i] - xmean) * (known_y[i] - ymean);
    den += Math.pow(known_x[i] - xmean, 2);
  }

  return num / den
}

/**
 * Returns the k-th smallest value in a data set.
 *
 * Category: Statistical
 *
 * @param {*} array An array or range of numerical data for which you want to determine the k-th smallest value.
 * @param {*} k The position (from the smallest) in the array or range of data to return.
 * @returns
 */
function SMALL(array, k) {
  array = parseNumberArray(flatten(array));
  k = parseNumber(k);

  if (anyIsError(array, k)) {
    return array
  }

  return array.sort((a, b) => a - b)[k - 1]
}

/**
 * Returns a normalized value.
 *
 * Category: Statistical
 *
 * @param {*} x The value you want to normalize.
 * @param {*} mean The arithmetic mean of the distribution.
 * @param {*} standard_dev The standard deviation of the distribution.
 * @returns
 */
function STANDARDIZE(x, mean, standard_dev) {
  x = parseNumber(x);
  mean = parseNumber(mean);
  standard_dev = parseNumber(standard_dev);

  if (anyIsError(x, mean, standard_dev)) {
    return value
  }

  return (x - mean) / standard_dev
}

const STDEV = {};

/**
 * Calculates standard deviation based on the entire population.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number arguments 2 to 254 corresponding to a population. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
STDEV.P = function () {
  const v = VAR.P.apply(this, arguments);
  let result = Math.sqrt(v);

  if (isNaN(result)) {
    result = num;
  }

  return result
};

/**
 * Estimates standard deviation based on a sample.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number arguments 2 to 254 corresponding to a sample of a population. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
STDEV.S = function () {
  const v = VAR.S.apply(this, arguments);
  const result = Math.sqrt(v);

  return result
};

/**
 * Estimates standard deviation based on a sample, including numbers, text, and logical values.
 *
 * Category: Statistical
 *
 * @param {*} args value1, value2, ... Value1 is required, subsequent values are optional. 1 to 255 values corresponding to a sample of a population. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
function STDEVA() {
  const v = VARA.apply(this, arguments);
  const result = Math.sqrt(v);

  return result
}

/**
 * Calculates standard deviation based on the entire population, including numbers, text, and logical values.
 *
 * Category: Statistical
 *
 * @param {*} args value1, value2, ... Value1 is required, subsequent values are optional. 1 to 255 values corresponding to a population. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
function STDEVPA() {
  const v = VARPA.apply(this, arguments);
  let result = Math.sqrt(v);

  if (isNaN(result)) {
    result = num;
  }

  return result
}

/**
 * Returns the standard error of the predicted y-value for each x in the regression.
 *
 * Category: Statistical
 *
 * @param {*} known_y An array or range of dependent data points.
 * @param {*} known_x An array or range of independent data points.
 * @returns
 */
function STEYX(known_y, known_x) {
  known_y = parseNumberArray(flatten(known_y));
  known_x = parseNumberArray(flatten(known_x));

  if (anyIsError(known_y, known_x)) {
    return value
  }

  const xmean = jstat.mean(known_x);
  const ymean = jstat.mean(known_y);
  const n = known_x.length;
  let lft = 0;
  let num = 0;
  let den = 0;

  for (let i = 0; i < n; i++) {
    lft += Math.pow(known_y[i] - ymean, 2);
    num += (known_x[i] - xmean) * (known_y[i] - ymean);
    den += Math.pow(known_x[i] - xmean, 2);
  }

  return Math.sqrt((lft - (num * num) / den) / (n - 2))
}

const T$1 = {};

/**
 * Returns the Percentage Points (probability) for the Student t-distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The numeric value at which to evaluate the distribution
 * @param {*} deg_freedom An integer indicating the number of degrees of freedom.
 * @param {*} cumulative A logical value that determines the form of the function. If cumulative is TRUE, T.DIST returns the cumulative distribution function; if FALSE, it returns the probability density function.
 * @returns
 */
T$1.DIST = (x, deg_freedom, cumulative) => {
  if (cumulative !== 1 && cumulative !== 2) {
    return num
  }

  return cumulative === 1 ? T$1.DIST.RT(x, deg_freedom) : T$1.DIST['2T'](x, deg_freedom)
};

/**
 * Returns the Percentage Points (probability) for the Student t-distribution
 *
 * Category: Statistical
 *
 * @param {*} x The numeric value at which to evaluate the distribution.
 * @param {*} deg_freedom An integer indicating the number of degrees of freedom.
 * @returns
 */
T$1.DIST['2T'] = function (x, deg_freedom) {
  if (arguments.length !== 2) {
    return na
  }

  if (x < 0 || deg_freedom < 1) {
    return num
  }

  if (typeof x !== 'number' || typeof deg_freedom !== 'number') {
    return value
  }

  return (1 - jstat.studentt.cdf(x, deg_freedom)) * 2
};

/**
 * Returns the Student's t-distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The numeric value at which to evaluate the distribution.
 * @param {*} deg_freedom An integer indicating the number of degrees of freedom.
 * @returns
 */
T$1.DIST.RT = function (x, deg_freedom) {
  if (arguments.length !== 2) {
    return na
  }

  if (x < 0 || deg_freedom < 1) {
    return num
  }

  if (typeof x !== 'number' || typeof deg_freedom !== 'number') {
    return value
  }

  return 1 - jstat.studentt.cdf(x, deg_freedom)
};

/**
 * Returns the t-value of the Student's t-distribution as a function of the probability and the degrees of freedom.
 *
 * Category: Statistical
 *
 * @param {*} probability The probability associated with the Student's t-distribution.
 * @param {*} deg_freedom The number of degrees of freedom with which to characterize the distribution.
 * @returns
 */
T$1.INV = (probability, deg_freedom) => {
  probability = parseNumber(probability);
  deg_freedom = parseNumber(deg_freedom);

  if (anyIsError(probability, deg_freedom)) {
    return value
  }

  return jstat.studentt.inv(probability, deg_freedom)
};

/**
 * Returns the inverse of the Student's t-distribution
 *
 * Category: Statistical
 *
 * @param {*} probability The probability associated with the Student's t-distribution.
 * @param {*} deg_freedom The number of degrees of freedom with which to characterize the distribution.
 * @returns
 */
T$1.INV['2T'] = (probability, deg_freedom) => {
  probability = parseNumber(probability);
  deg_freedom = parseNumber(deg_freedom);

  if (probability <= 0 || probability > 1 || deg_freedom < 1) {
    return num
  }

  if (anyIsError(probability, deg_freedom)) {
    return value
  }

  return Math.abs(jstat.studentt.inv(probability / 2, deg_freedom))
};

// The algorithm can be found here:
// http://www.chem.uoa.gr/applets/AppletTtest/Appl_Ttest2.html
/**
 * Returns the probability associated with a Student's t-test.
 *
 * Category: Statistical
 *
 * @param {*} array1 The first data set.
 * @param {*} array2 The second data set.
 * @returns
 */
T$1.TEST = (array1, array2) => {
  array1 = parseNumberArray(flatten(array1));
  array2 = parseNumberArray(flatten(array2));

  if (anyIsError(array1, array2)) {
    return value
  }

  const mean_x = jstat.mean(array1);
  const mean_y = jstat.mean(array2);
  let s_x = 0;
  let s_y = 0;
  let i;

  for (i = 0; i < array1.length; i++) {
    s_x += Math.pow(array1[i] - mean_x, 2);
  }

  for (i = 0; i < array2.length; i++) {
    s_y += Math.pow(array2[i] - mean_y, 2);
  }

  s_x = s_x / (array1.length - 1);
  s_y = s_y / (array2.length - 1);

  const t = Math.abs(mean_x - mean_y) / Math.sqrt(s_x / array1.length + s_y / array2.length);

  return T$1.DIST['2T'](t, array1.length + array2.length - 2)
};

/**
 * Returns values along a linear trend.
 *
 * Category: Statistical
 *
 * @param {*} known_ys The set of y-values you already know in the relationship y = mx + b
 * @param {*} known_xs An optional set of x-values that you may already know in the relationship y = mx + b
 * @param {*} new_xs Optional. New x-values for which you want TREND to return corresponding y-values.
 * @returns
 */
function TREND(known_ys, known_xs, new_xs) {
  known_ys = parseNumberArray(flatten(known_ys));
  known_xs = parseNumberArray(flatten(known_xs));
  new_xs = parseNumberArray(flatten(new_xs));

  if (anyIsError(known_ys, known_xs, new_xs)) {
    return value
  }

  const linest = LINEST(known_ys, known_xs);
  const m = linest[0];
  const b = linest[1];
  const result = [];

  new_xs.forEach((x) => {
    result.push(m * x + b);
  });

  return result
}

/**
 * Returns the mean of the interior of a data set.
 *
 * Category: Statistical
 *
 * @param {*} array The array or range of values to trim and average.
 * @param {*} percent The fractional number of data points to exclude from the calculation. For example, if percent = 0.2, 4 points are trimmed from a data set of 20 points (20 x 0.2): 2 from the top and 2 from the bottom of the set.
 * @returns
 */
function TRIMMEAN(range, percent) {
  range = parseNumberArray(flatten(range));
  percent = parseNumber(percent);

  if (anyIsError(range, percent)) {
    return value
  }

  const trim = FLOOR(range.length * percent, 2) / 2;

  return jstat.mean(
    initial(
      rest(
        range.sort((a, b) => a - b),
        trim
      ),
      trim
    )
  )
}

const VAR = {};

/**
 * Calculates variance based on the entire population.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number arguments 2 to 254 corresponding to a population.
 * @returns
 */
VAR.P = function () {
  const range = numbers(flatten(arguments));
  const n = range.length;
  let sigma = 0;
  const mean = AVERAGE(range);
  let result;

  for (let i = 0; i < n; i++) {
    sigma += Math.pow(range[i] - mean, 2);
  }

  result = sigma / n;

  if (isNaN(result)) {
    result = num;
  }

  return result
};

/**
 * Estimates variance based on a sample.
 *
 * Category: Statistical
 *
 * @param {*} args number1, number2, ... Number arguments 2 to 254 corresponding to a sample of a population.
 * @returns
 */
VAR.S = function () {
  const range = numbers(flatten(arguments));
  const n = range.length;
  let sigma = 0;
  const mean = AVERAGE(range);

  for (let i = 0; i < n; i++) {
    sigma += Math.pow(range[i] - mean, 2);
  }

  return sigma / (n - 1)
};

/**
 * Estimates variance based on a sample, including numbers, text, and logical values.
 *
 * Category: Statistical
 *
 * @param {*} args value1, value2, ... Value1 is required, subsequent values are optional. 1 to 255 value arguments corresponding to a sample of a population.
 * @returns
 */
function VARA() {
  const range = flatten(arguments);
  const n = range.length;
  let sigma = 0;
  let count = 0;
  const mean = AVERAGEA(range);

  for (let i = 0; i < n; i++) {
    const el = range[i];

    if (typeof el === 'number') {
      sigma += Math.pow(el - mean, 2);
    } else if (el === true) {
      sigma += Math.pow(1 - mean, 2);
    } else {
      sigma += Math.pow(0 - mean, 2);
    }

    if (el !== null) {
      count++;
    }
  }

  return sigma / (count - 1)
}

/**
 * Calculates variance based on the entire population, including numbers, text, and logical values.
 *
 * Category: Statistical
 *
 * @param {*} args value1, value2, ... Value1 is required, subsequent values are optional. 1 to 255 value arguments corresponding to a population.
 * @returns
 */
function VARPA() {
  const range = flatten(arguments);
  const n = range.length;
  let sigma = 0;
  let count = 0;
  const mean = AVERAGEA(range);
  let result;

  for (let i = 0; i < n; i++) {
    const el = range[i];

    if (typeof el === 'number') {
      sigma += Math.pow(el - mean, 2);
    } else if (el === true) {
      sigma += Math.pow(1 - mean, 2);
    } else {
      sigma += Math.pow(0 - mean, 2);
    }

    if (el !== null) {
      count++;
    }
  }

  result = sigma / count;

  if (isNaN(result)) {
    result = num;
  }

  return result
}

const WEIBULL = {};

/**
 * Returns the Weibull distribution.
 *
 * Category: Statistical
 *
 * @param {*} x The value at which to evaluate the function.
 * @param {*} alpha A parameter to the distribution.
 * @param {*} beta A parameter to the distribution.
 * @param {*} cumulative Determines the form of the function.
 * @returns
 */
WEIBULL.DIST = (x, alpha, beta, cumulative) => {
  x = parseNumber(x);
  alpha = parseNumber(alpha);
  beta = parseNumber(beta);

  if (anyIsError(x, alpha, beta)) {
    return value
  }

  return cumulative
    ? 1 - Math.exp(-Math.pow(x / beta, alpha))
    : (Math.pow(x, alpha - 1) * Math.exp(-Math.pow(x / beta, alpha)) * alpha) / Math.pow(beta, alpha)
};

const Z = {};

/**
 * Returns the one-tailed probability-value of a z-test.
 *
 * Category: Statistical
 *
 * @param {*} array The array or range of data against which to test x.
 * @param {*} x The value to test.
 * @param {*} sigma Optional. The population (known) standard deviation. If omitted, the sample standard deviation is used.
 * @returns
 */
Z.TEST = (array, x, sigma) => {
  array = parseNumberArray(flatten(array));
  x = parseNumber(x);

  if (anyIsError(array, x)) {
    return value
  }

  sigma = sigma || STDEV.S(array);
  const n = array.length;

  return 1 - NORM.S.DIST((AVERAGE(array) - x) / (sigma / Math.sqrt(n)), true)
};

/**
 * Returns the absolute value of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The real number of which you want the absolute value.
 * @returns
 */
function ABS(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  const result = Math.abs(number);

  return result
}

/**
 * Returns the arccosine of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The cosine of the angle you want and must be from -1 to 1.
 * @returns
 */
function ACOS(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  let result = Math.acos(number);

  if (isNaN(result)) {
    result = num;
  }

  return result
}

/**
 * Returns the inverse hyperbolic cosine of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Any real number equal to or greater than 1.
 * @returns
 */
function ACOSH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  let result = Math.log(number + Math.sqrt(number * number - 1));

  if (isNaN(result)) {
    result = num;
  }

  return result
}

/**
 * Returns the arccotangent of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Number is the cotangent of the angle you want. This must be a real number.
 * @returns
 */
function ACOT(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  const result = Math.atan(1 / number);

  return result
}

/**
 * Returns the hyperbolic arccotangent of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The absolute value of Number must be greater than 1.
 * @returns
 */
function ACOTH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  let result = 0.5 * Math.log((number + 1) / (number - 1));

  if (isNaN(result)) {
    result = num;
  }

  return result
}

// TODO: use options
/**
 * Returns an aggregate in a list or database.
 *
 * Category: Math and trigonometry
 *
 * @param {*} function_num A number 1 to 19 that specifies which function to use.
 * @param {*} options A numerical value that determines which values to ignore in the evaluation range for the function. Note: The function will not ignore hidden rows, nested subtotals or nested aggregates if the array argument includes a calculation, for example: =AGGREGATE(14,3,A1:A100*(A1:A100>0),1)
 * @param {*} ref1 The first numeric argument for functions that take multiple numeric arguments for which you want the aggregate value.
 * @param {*} ref2 Optional. Numeric arguments 2 to 253 for which you want the aggregate value. For functions that take an array, ref1 is an array, an array formula, or a reference to a range of values for which you want the aggregate value. Ref2 is a second argument that is required for certain functions.
 * @returns
 */
function AGGREGATE(function_num, options, ref1, ref2) {
  function_num = parseNumber(function_num);
  options = parseNumber(function_num);

  if (anyIsError(function_num, options)) {
    return value
  }

  switch (function_num) {
    case 1:
      return AVERAGE(ref1)
    case 2:
      return COUNT(ref1)
    case 3:
      return COUNTA(ref1)
    case 4:
      return MAX(ref1)
    case 5:
      return MIN(ref1)
    case 6:
      return PRODUCT(ref1)
    case 7:
      return STDEV.S(ref1)
    case 8:
      return STDEV.P(ref1)
    case 9:
      return SUM(ref1)
    case 10:
      return VAR.S(ref1)
    case 11:
      return VAR.P(ref1)
    case 12:
      return MEDIAN(ref1)
    case 13:
      return MODE.SNGL(ref1)
    case 14:
      return LARGE(ref1, ref2)
    case 15:
      return SMALL(ref1, ref2)
    case 16:
      return PERCENTILE.INC(ref1, ref2)
    case 17:
      return QUARTILE.INC(ref1, ref2)
    case 18:
      return PERCENTILE.EXC(ref1, ref2)
    case 19:
      return QUARTILE.EXC(ref1, ref2)
  }
}

/**
 * Converts a Roman number to Arabic, as a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} text A string enclosed in quotation marks, an empty string (""), or a reference to a value containing text.
 * @returns
 */
function ARABIC(text) {
  if (text === undefined || text === null) {
    return 0
  }

  if (text instanceof Error) {
    return text
  }

  // Credits: Rafa? Kukawski
  if (!/^M*(?:D?C{0,3}|C[MD])(?:L?X{0,3}|X[CL])(?:V?I{0,3}|I[XV])$/.test(text)) {
    return value
  }

  let r = 0;
  text.replace(/[MDLV]|C[MD]?|X[CL]?|I[XV]?/g, (i) => {
    r += {
      M: 1000,
      CM: 900,
      D: 500,
      CD: 400,
      C: 100,
      XC: 90,
      L: 50,
      XL: 40,
      X: 10,
      IX: 9,
      V: 5,
      IV: 4,
      I: 1
    }[i];
  });

  return r
}

/**
 * Returns the arcsine of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The sine of the angle you want and must be from -1 to 1.
 * @returns
 */
function ASIN(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  let result = Math.asin(number);

  if (isNaN(result)) {
    result = num;
  }

  return result
}

/**
 * Returns the inverse hyperbolic sine of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Any real number.
 * @returns
 */
function ASINH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return Math.log(number + Math.sqrt(number * number + 1))
}

/**
 * Returns the arctangent of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The tangent of the angle you want.
 * @returns
 */
function ATAN(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return Math.atan(number)
}

/**
 * Returns the arctangent from x- and y-coordinates.
 *
 * Category: Math and trigonometry
 *
 * @param {*} x_num The x-coordinate of the point.
 * @param {*} y_num The y-coordinate of the point.
 * @returns
 */
function ATAN2(x_num, y_num) {
  x_num = parseNumber(x_num);
  y_num = parseNumber(y_num);
  const anyError$1 = anyError(x_num, y_num);

  if (anyError$1) {
    return anyError$1
  }

  return Math.atan2(x_num, y_num)
}

/**
 * Returns the inverse hyperbolic tangent of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Any real number between 1 and -1.
 * @returns
 */
function ATANH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  let result = Math.log((1 + number) / (1 - number)) / 2;

  if (isNaN(result)) {
    result = num;
  }

  return result
}

/**
 * Converts a number into a text representation with the given radix (base).
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The number that you want to convert. Must be an integer greater than or equal to 0 and less than 2^53.
 * @param {*} radix The base radix that you want to convert the number into. Must be an integer greater than or equal to 2 and less than or equal to 36.
 * @param {*} min_length Optional. The minimum length of the returned string. Must be an integer greater than or equal to 0.
 * @returns
 */
function BASE(number, radix, min_length) {
  number = parseNumber(number);
  radix = parseNumber(radix);
  min_length = parseNumber(min_length);
  const anyError$1 = anyError(number, radix, min_length);

  if (anyError$1) {
    return anyError$1
  }

  if (radix === 0) {
    return num
  }

  const result = number.toString(radix);

  return new Array(Math.max(min_length + 1 - result.length, 0)).join('0') + result
}

/**
 * Rounds a number to the nearest integer or to the nearest multiple of significance.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The value you want to round.
 * @param {*} significance The multiple to which you want to round.
 * @param {*} mode Optional. For negative numbers, controls whether Number is rounded toward or away from zero.
 * @returns
 */
function CEILING(number, significance, mode) {
  number = parseNumber(number);
  significance = parseNumber(significance);
  mode = parseNumber(mode);
  const anyError$1 = anyError(number, significance, mode);

  if (anyError$1) {
    return anyError$1
  }

  if (significance === 0) {
    return 0
  }

  significance = Math.abs(significance);
  const precision = -Math.floor(Math.log(significance) / Math.log(10));

  if (number >= 0) {
    return ROUND(Math.ceil(number / significance) * significance, precision)
  } else {
    if (mode === 0) {
      return -ROUND(Math.floor(Math.abs(number) / significance) * significance, precision)
    } else {
      return -ROUND(Math.ceil(Math.abs(number) / significance) * significance, precision)
    }
  }
}

CEILING.MATH = CEILING;

CEILING.PRECISE = CEILING;

/**
 * Returns the number of combinations for a given number of objects.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The number of items.
 * @param {*} number_chosen The number of items in each combination.
 * @returns
 */
function COMBIN(number, number_chosen) {
  number = parseNumber(number);
  number_chosen = parseNumber(number_chosen);
  const anyError$1 = anyError(number, number_chosen);

  if (anyError$1) {
    return anyError$1
  }

  if (number < number_chosen) {
    return num
  }

  return FACT(number) / (FACT(number_chosen) * FACT(number - number_chosen))
}

/**
 * Returns the number of combinations with repetitions for a given number of items.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Must be greater than or equal to 0, and greater than or equal to Number_chosen. Non-integer values are truncated.
 * @param {*} number_chosen Must be greater than or equal to 0. Non-integer values are truncated.
 * @returns
 */
function COMBINA(number, number_chosen) {
  number = parseNumber(number);
  number_chosen = parseNumber(number_chosen);
  const anyError$1 = anyError(number, number_chosen);

  if (anyError$1) {
    return anyError$1
  }

  if (number < number_chosen) {
    return num
  }

  return number === 0 && number_chosen === 0 ? 1 : COMBIN(number + number_chosen - 1, number - 1)
}

/**
 * Returns the cosine of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The angle in radians for which you want the cosine.
 * @returns
 */
function COS(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return Math.cos(number)
}

/**
 * Returns the hyperbolic cosine of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Any real number for which you want to find the hyperbolic cosine.
 * @returns
 */
function COSH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return (Math.exp(number) + Math.exp(-number)) / 2
}

/**
 * Returns the hyperbolic cosine of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The angle in radians for which you want the cotangent.
 * @returns
 */
function COT(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  if (number === 0) {
    return div0
  }

  return 1 / Math.tan(number)
}

/**
 * Returns the cotangent of an angle.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number
 * @returns
 */
function COTH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  if (number === 0) {
    return div0
  }

  const e2 = Math.exp(2 * number);

  return (e2 + 1) / (e2 - 1)
}

/**
 * Returns the cosecant of an angle.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number
 * @returns
 */
function CSC(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  if (number === 0) {
    return div0
  }

  return 1 / Math.sin(number)
}

/**
 * Returns the hyperbolic cosecant of an angle.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number
 * @returns
 */
function CSCH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  if (number === 0) {
    return div0
  }

  return 2 / (Math.exp(number) - Math.exp(-number))
}

/**
 * Converts a text representation of a number in a given base into a decimal number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} text
 * @param {*} radix Radix must be an integer.
 * @returns
 */
function DECIMAL(text, radix) {
  if (arguments.length < 1) {
    return value
  }

  text = parseNumber(text);
  radix = parseNumber(radix);
  const anyError$1 = anyError(text, radix);

  if (anyError$1) {
    return anyError$1
  }

  if (radix === 0) {
    return num
  }

  return parseInt(text, radix)
}

/**
 * Converts radians to degrees.
 *
 * Category: Math and trigonometry
 *
 * @param {*} angle The angle in radians that you want to convert.
 * @returns
 */
function DEGREES(angle) {
  angle = parseNumber(angle);

  if (angle instanceof Error) {
    return angle
  }

  return (angle * 180) / Math.PI
}

/**
 * Rounds a number up to the nearest even integer.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The value to round.
 * @returns
 */
function EVEN(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return CEILING(number, -2, -1)
}

/**
 * Returns e raised to the power of a given number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The exponent applied to the base e.
 * @returns
 */
function EXP(number) {
  if (arguments.length < 1) {
    return na
  }

  if (arguments.length > 1) {
    return error
  }

  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  number = Math.exp(number);

  return number
}

const MEMOIZED_FACT = [];
/**
 * Returns the factorial of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The nonnegative number for which you want the factorial. If number is not an integer, it is truncated.
 * @returns
 */
function FACT(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  const n = Math.floor(number);

  if (n === 0 || n === 1) {
    return 1
  } else if (MEMOIZED_FACT[n] > 0) {
    return MEMOIZED_FACT[n]
  } else {
    MEMOIZED_FACT[n] = FACT(n - 1) * n;

    return MEMOIZED_FACT[n]
  }
}

/**
 * Returns the double factorial of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The value for which to return the double factorial. If number is not an integer, it is truncated.
 * @returns
 */
function FACTDOUBLE(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  const n = Math.floor(number);

  return n <= 0 ? 1 : n * FACTDOUBLE(n - 2)
}

/**
 * Rounds a number down, toward zero.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The numeric value you want to round.
 * @param {*} significance The multiple to which you want to round.
 * @returns
 */
function FLOOR(number, significance) {
  number = parseNumber(number);
  significance = parseNumber(significance);
  const anyError$1 = anyError(number, significance);

  if (anyError$1) {
    return anyError$1
  }

  if (significance === 0) {
    return 0
  }

  if (!(number >= 0 && significance > 0) && !(number <= 0 && significance < 0)) {
    return num
  }

  significance = Math.abs(significance);
  const precision = -Math.floor(Math.log(significance) / Math.log(10));

  return number >= 0
    ? ROUND(Math.floor(number / significance) * significance, precision)
    : -ROUND(Math.ceil(Math.abs(number) / significance), precision)
}

// TODO: Verify

/**
 * Rounds a number down, to the nearest integer or to the nearest multiple of significance.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The number to be rounded down.
 * @param {*} significance Optional. The multiple to which you want to round.
 * @param {*} mode Optional. The direction (toward or away from 0) to round negative numbers.
 * @returns
 */
FLOOR.MATH = (number, significance, mode) => {
  if (significance instanceof Error) {
    return significance
  }

  significance = significance === undefined ? 0 : significance;

  number = parseNumber(number);
  significance = parseNumber(significance);
  mode = parseNumber(mode);
  const anyError$1 = anyError(number, significance, mode);

  if (anyError$1) {
    return anyError$1
  }

  if (significance === 0) {
    return 0
  }

  significance = significance ? Math.abs(significance) : 1;
  const precision = -Math.floor(Math.log(significance) / Math.log(10));

  if (number >= 0) {
    return ROUND(Math.floor(number / significance) * significance, precision)
  } else if (mode === 0 || mode === undefined) {
    return -ROUND(Math.ceil(Math.abs(number) / significance) * significance, precision)
  }

  return -ROUND(Math.floor(Math.abs(number) / significance) * significance, precision)
};

// Deprecated

/**
 * Rounds a number the nearest integer or to the nearest multiple of significance. Regardless of the sign of the number, the number is rounded up.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The value to be rounded.
 * @param {*} significance Optional. The multiple to which number is to be rounded. If significance is omitted, its default value is 1.
 * @returns
 */
FLOOR.PRECISE = FLOOR['MATH'];

// adapted http://rosettacode.org/wiki/Greatest_common_divisor#JavaScript
/**
 * Returns the greatest common divisor.
 *
 * Category: Math and trigonometry
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 values. If any value is not an integer, it is truncated.
 * @returns
 */
function GCD() {
  const range = parseNumberArray(flatten(arguments));

  if (range instanceof Error) {
    return range
  }

  const n = range.length;
  const r0 = range[0];
  let x = r0 < 0 ? -r0 : r0;

  for (let i = 1; i < n; i++) {
    const ri = range[i];
    let y = ri < 0 ? -ri : ri;

    while (x && y) {
      if (x > y) {
        x %= y;
      } else {
        y %= x;
      }
    }

    x += y;
  }

  return x
}

/**
 * Rounds a number down to the nearest integer.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The real number you want to round down to an integer.
 * @returns
 */
function INT(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return Math.floor(number)
}

// TODO: verify
const ISO = {
  CEILING: CEILING
};

/**
 * Returns the least common multiple.
 *
 * Category: Math and trigonometry
 *
 * @param {*} args number1, number2,... Number1 is required, subsequent numbers are optional. 1 to 255 values for which you want the least common multiple. If value is not an integer, it is truncated.
 * @returns
 */
function LCM() {
  // Credits: Jonas Raoni Soares Silva
  const o = parseNumberArray(flatten(arguments));

  if (o instanceof Error) {
    return o
  }

  for (var i, j, n, d, r = 1; (n = o.pop()) !== undefined; ) {
    if (n === 0) {
      return 0
    }

    while (n > 1) {
      if (n % 2) {
        for (i = 3, j = Math.floor(Math.sqrt(n)); i <= j && n % i; i += 2) {
          // empty
        }

        d = i <= j ? i : n;
      } else {
        d = 2;
      }

      for (n /= d, r *= d, i = o.length; i; o[--i] % d === 0 && (o[i] /= d) === 1 && o.splice(i, 1)) {
        // empty
      }
    }
  }

  return r
}

/**
 * Returns the natural logarithm of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The positive real number for which you want the natural logarithm.
 * @returns
 */
function LN(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  if (number === 0) {
    return num
  }

  return Math.log(number)
}

/**
 * Returns the logarithm of a number to a specified base.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The positive real number for which you want the logarithm.
 * @param {*} base Optional. The base of the logarithm. If base is omitted, it is assumed to be 10.
 * @returns
 */
function LOG(number, base) {
  number = parseNumber(number);
  base = parseNumber(base);
  const anyError$1 = anyError(number, base);

  if (anyError$1) {
    return anyError$1
  }

  if (number === 0 || base === 0) {
    return num
  }

  return Math.log(number) / Math.log(base)
}

/**
 * Returns the base-10 logarithm of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The positive real number for which you want the base-10 logarithm.
 * @returns
 */
function LOG10(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  if (number === 0) {
    return num
  }

  return Math.log(number) / Math.log(10)
}

/**
 * Returns the matrix product of two arrays. The result is an array with the same number of rows as array1 and the same number of columns as array2.
 *
 * Category: Math and trigonometry
 *
 * @param {*} array1 Required. 1st array you want to multiply.
 * @param {*} array2 Required. 2nd array you want to multiply.
 * @returns
 */
function MMULT(array1, array2) {
  if (
    //Arguments are not arrays
    !Array.isArray(array1) ||
    !Array.isArray(array2) ||
    // There are empty arrays
    array1.some((el) => !el.length) ||
    array2.some((el) => !el.length) ||
    // Not all array elements are numbers
    flattenShallow(array1).some((el) => typeof el !== 'number') ||
    flattenShallow(array2).some((el) => typeof el !== 'number') ||
    // Number of columns in array1 is different from the number of rows in array2
    array1[0].length !== array2.length
  ) {
    return value
  }

  const matrix = Array(array1.length)
    .fill(0)
    .map(() => Array(array2[0].length).fill(0));

  return matrix.map((row, i) => row.map((_, j) => array1[i].reduce((sum, el, k) => sum + el * array2[k][j], 0)))
}

/**
 * Returns the remainder from division.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The number for which you want to find the remainder.
 * @param {*} divisor The number by which you want to divide number.
 * @returns
 */
function MOD(number, divisor) {
  number = parseNumber(number);
  divisor = parseNumber(divisor);
  const anyError$1 = anyError(number, divisor);

  if (anyError$1) {
    return anyError$1
  }

  if (divisor === 0) {
    return div0
  }

  let modulus = Math.abs(number % divisor);
  modulus = number < 0 ? divisor - modulus : modulus;

  return divisor > 0 ? modulus : -modulus
}

/**
 * Returns a number rounded to the desired multiple.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The value to round.
 * @param {*} multiple The multiple to which you want to round number.
 * @returns
 */
function MROUND(number, multiple) {
  number = parseNumber(number);
  multiple = parseNumber(multiple);
  const anyError$1 = anyError(number, multiple);

  if (anyError$1) {
    return anyError$1
  }

  if (number * multiple === 0) {
    return 0
  }

  if (number * multiple < 0) {
    return num
  }

  return Math.round(number / multiple) * multiple
}

/**
 * Returns the multinomial of a set of numbers.
 *
 * Category: Math and trigonometry
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 values for which you want the multinomial.
 * @returns
 */
function MULTINOMIAL() {
  const args = parseNumberArray(flatten(arguments));

  if (args instanceof Error) {
    return args
  }

  let sum = 0;
  let divisor = 1;

  for (let i = 0; i < args.length; i++) {
    sum += args[i];
    divisor *= FACT(args[i]);
  }

  return FACT(sum) / divisor
}

/**
 * Returns the unit matrix for the specified dimension.
 *
 * Category: Math and trigonometry
 *
 * @param {*} dimension Required. Dimension is an integer specifying the dimension of the unit matrix that you want to return. It returns an array. The dimension has to be greater than zero.
 * @returns
 */
function MUNIT(dimension) {
  if (arguments.length > 1) {
    return na
  }

  dimension = parseInt(dimension);

  if (!dimension || dimension <= 0) {
    return value
  }

  return Array(dimension)
    .fill(0)
    .map(() => Array(dimension).fill(0))
    .map((el, i) => {
      el[i] = 1;
      return el
    })
}

/**
 * Rounds a number up to the nearest odd integer.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number: The value to round.
 * @returns
 */
function ODD(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  let temp = Math.ceil(Math.abs(number));
  temp = temp & 1 ? temp : temp + 1;

  return number >= 0 ? temp : -temp
}

/**
 * Returns the value of pi.
 *
 * Category: Math and trigonometry
 *
 * @returns
 */
function PI() {
  return Math.PI
}

/**
 * Returns the result of a number raised to a power.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The base number. It can be any real number.
 * @param {*} power The exponent to which the base number is raised.
 * @returns
 */
function POWER(number, power) {
  number = parseNumber(number);
  power = parseNumber(power);
  const anyError$1 = anyError(number, power);

  if (anyError$1) {
    return anyError$1
  }

  if (number === 0 && power === 0) {
    return num
  }

  const result = Math.pow(number, power);

  if (isNaN(result)) {
    return num
  }

  return result
}

/**
 * Multiplies its arguments.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number1 The first number or range that you want to multiply.
 * @param {*} args number2, ... Optional. Additional numbers or ranges that you want to multiply, up to a maximum of 255 arguments.
 * @returns
 */
function PRODUCT() {
  const flatArguments = flatten(arguments);
  const flatArgumentsDefined = flatArguments.filter((arg) => arg !== undefined && arg !== null);

  if (flatArgumentsDefined.length === 0) {
    return 0
  }

  const args = parseNumberArray(flatArgumentsDefined);

  if (args instanceof Error) {
    return args
  }

  let result = 1;

  for (let i = 0; i < args.length; i++) {
    result *= args[i];
  }

  return result
}

/**
 * Returns the integer portion of a division.
 *
 * Category: Math and trigonometry
 *
 * @param {*} numerator The dividend.
 * @param {*} denominator The divisor.
 * @returns
 */
function QUOTIENT(numerator, denominator) {
  numerator = parseNumber(numerator);
  denominator = parseNumber(denominator);
  const anyError$1 = anyError(numerator, denominator);

  if (anyError$1) {
    return anyError$1
  }

  return parseInt(numerator / denominator, 10)
}

/**
 * Converts degrees to radians.
 *
 * Category: Math and trigonometry
 *
 * @param {*} angle An angle in degrees that you want to convert.
 * @returns
 */
function RADIANS(angle) {
  angle = parseNumber(angle);

  if (angle instanceof Error) {
    return angle
  }

  return (angle * Math.PI) / 180
}

/**
 * Returns a random number between 0 and 1.
 *
 * Category: Math and trigonometry
 *
 * @returns
 */
function RAND() {
  return Math.random()
}

/**
 * Returns a random number between the numbers you specify.
 *
 * Category: Math and trigonometry
 *
 * @param {*} bottom The smallest integer RANDBETWEEN will return.
 * @param {*} top The largest integer RANDBETWEEN will return.
 * @returns
 */
function RANDBETWEEN(bottom, top) {
  bottom = parseNumber(bottom);
  top = parseNumber(top);
  const anyError$1 = anyError(bottom, top);

  if (anyError$1) {
    return anyError$1
  }
  // Creative Commons Attribution 3.0 License
  // Copyright (c) 2012 eqcode

  return bottom + Math.ceil((top - bottom + 1) * Math.random()) - 1
}

// TODO
/**
 * Converts an arabic numeral to roman, as text.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The Arabic numeral you want converted.
 * @returns
 */
function ROMAN(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  // The MIT License
  // Copyright (c) 2008 Steven Levithan
  const digits = String(number).split('');
  const key = [
    '',
    'C',
    'CC',
    'CCC',
    'CD',
    'D',
    'DC',
    'DCC',
    'DCCC',
    'CM',
    '',
    'X',
    'XX',
    'XXX',
    'XL',
    'L',
    'LX',
    'LXX',
    'LXXX',
    'XC',
    '',
    'I',
    'II',
    'III',
    'IV',
    'V',
    'VI',
    'VII',
    'VIII',
    'IX'
  ];
  let roman = '';
  let i = 3;

  while (i--) {
    roman = (key[+digits.pop() + i * 10] || '') + roman;
  }

  return new Array(+digits.join('') + 1).join('M') + roman
}

/**
 * Rounds a number to a specified number of digits.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The number that you want to round.
 * @param {*} num_digits The number of digits to which you want to round the number argument.
 * @returns
 */
function ROUND(number, num_digits) {
  number = parseNumber(number);
  num_digits = parseNumber(num_digits);
  const anyError$1 = anyError(number, num_digits);

  if (anyError$1) {
    return anyError$1
  }

  return Number(Math.round(Number(number + 'e' + num_digits)) + 'e' + num_digits * -1)
}

/**
 * Rounds a number down, toward zero.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Any real number that you want rounded down.
 * @param {*} num_digits The number of digits to which you want to round number.
 * @returns
 */
function ROUNDDOWN(number, num_digits) {
  number = parseNumber(number);
  num_digits = parseNumber(num_digits);
  const anyError$1 = anyError(number, num_digits);

  if (anyError$1) {
    return anyError$1
  }

  const sign = number > 0 ? 1 : -1;

  return (sign * Math.floor(Math.abs(number) * Math.pow(10, num_digits))) / Math.pow(10, num_digits)
}

/**
 * Rounds a number up, away from zero.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Any real number that you want rounded up.
 * @param {*} num_digits The number of digits to which you want to round number.
 * @returns
 */
function ROUNDUP(number, num_digits) {
  number = parseNumber(number);
  num_digits = parseNumber(num_digits);
  const anyError$1 = anyError(number, num_digits);

  if (anyError$1) {
    return anyError$1
  }

  const sign = number > 0 ? 1 : -1;

  return (sign * Math.ceil(Math.abs(number) * Math.pow(10, num_digits))) / Math.pow(10, num_digits)
}

/**
 * Returns the secant of an angle.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The angle in radians for which you want the secant.
 * @returns
 */
function SEC(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return 1 / Math.cos(number)
}

/**
 * Returns the hyperbolic secant of an angle.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The angle in radians for which you want the hyperbolic secant.
 * @returns
 */
function SECH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return 2 / (Math.exp(number) + Math.exp(-number))
}

/**
 * Returns the sum of a power series based on the formula.
 *
 * Category: Math and trigonometry
 *
 * @param {*} x The input value to the power series.
 * @param {*} n The initial power to which you want to raise x.
 * @param {*} m The step by which to increase n for each term in the series.
 * @param {*} coefficients A set of coefficients by which each successive power of x is multiplied. The number of values in coefficients determines the number of terms in the power series. For example, if there are three values in coefficients, then there will be three terms in the power series.
 * @returns
 */
function SERIESSUM(x, n, m, coefficients) {
  x = parseNumber(x);
  n = parseNumber(n);
  m = parseNumber(m);
  coefficients = parseNumberArray(coefficients);

  if (anyIsError(x, n, m, coefficients)) {
    return value
  }

  let result = coefficients[0] * Math.pow(x, n);

  for (let i = 1; i < coefficients.length; i++) {
    result += coefficients[i] * Math.pow(x, n + i * m);
  }

  return result
}

/**
 * Returns the sign of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Any real number.
 * @returns
 */
function SIGN(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  if (number < 0) {
    return -1
  } else if (number === 0) {
    return 0
  } else {
    return 1
  }
}

/**
 * Returns the sine of the given angle.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The angle in radians for which you want the sine.
 * @returns
 */
function SIN(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return Math.sin(number)
}

/**
 * Returns the hyperbolic sine of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Any real number.
 * @returns
 */
function SINH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return (Math.exp(number) - Math.exp(-number)) / 2
}

/**
 * Returns a positive square root.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The number for which you want the square root.
 * @returns
 */
function SQRT(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  if (number < 0) {
    return num
  }

  return Math.sqrt(number)
}

/**
 * Returns the square root of (number * pi).
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The number by which pi is multiplied.
 * @returns
 */
function SQRTPI(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return Math.sqrt(number * Math.PI)
}

/**
 * Returns a subtotal in a list or database.
 *
 * Category: Math and trigonometry
 *
 * @param {*} function_num The number 1-11 or 101-111 that specifies the function to use for the subtotal. 1-11 includes manually-hidden rows, while 101-111 excludes them; filtered-out values are always excluded.
 * @param {*} ref1 The first named range or reference for which you want the subtotal.
 * @returns
 */
function SUBTOTAL(function_num, ref1) {
  function_num = parseNumber(function_num);

  if (function_num instanceof Error) {
    return function_num
  }

  switch (function_num) {
    case 1:
      return AVERAGE(ref1)
    case 2:
      return COUNT(ref1)
    case 3:
      return COUNTA(ref1)
    case 4:
      return MAX(ref1)
    case 5:
      return MIN(ref1)
    case 6:
      return PRODUCT(ref1)
    case 7:
      return STDEV.S(ref1)
    case 8:
      return STDEV.P(ref1)
    case 9:
      return SUM(ref1)
    case 10:
      return VAR.S(ref1)
    case 11:
      return VAR.P(ref1)
    // no hidden values for us
    case 101:
      return AVERAGE(ref1)
    case 102:
      return COUNT(ref1)
    case 103:
      return COUNTA(ref1)
    case 104:
      return MAX(ref1)
    case 105:
      return MIN(ref1)
    case 106:
      return PRODUCT(ref1)
    case 107:
      return STDEV.S(ref1)
    case 108:
      return STDEV.P(ref1)
    case 109:
      return SUM(ref1)
    case 110:
      return VAR.S(ref1)
    case 111:
      return VAR.P(ref1)
  }
}

/**
 * Adds its arguments.
 *
 * Category: Math and trigonometry
 *
 * @returns
 */
function SUM() {
  let result = 0;

  arrayEach(argsToArray(arguments), (value) => {
    if (result instanceof Error) {
      return false
    } else if (value instanceof Error) {
      result = value;
    } else if (typeof value === 'number') {
      result += value;
    } else if (typeof value === 'string') {
      const parsed = parseFloat(value);

      !isNaN(parsed) && (result += parsed);
    } else if (Array.isArray(value)) {
      const inner_result = SUM.apply(null, value);

      if (inner_result instanceof Error) {
        result = inner_result;
      } else {
        result += inner_result;
      }
    }
  });

  return result
}

/**
 * Adds the values specified by a given criteria.
 *
 * Category: Math and trigonometry
 *
 * @param {*} range The range of values that you want evaluated by criteria. Cells in each range must be numbers or names, arrays, or references that contain numbers. Blank and text values are ignored.
 * @param {*} criteria The criteria in the form of a number, expression, a value reference, text, or a function that defines which values will be added.
 * @param {*} sum_range Optional. The actual values to add, if you want to add values other than those specified in the range argument. If the sum_range argument is omitted, Excel adds the values that are specified in the range argument (the same values to which the criteria is applied). Sum_range should be the same size and shape as range. If it isn't, performance may suffer, and the formula will sum a range of values that starts with the first value in sum_range but has the same dimensions as range.
 * @returns
 */
function SUMIF(range, criteria, sum_range) {
  range = flatten(range);

  sum_range = sum_range ? flatten(sum_range) : range;

  if (range instanceof Error) {
    return range
  }

  if (criteria === undefined || criteria === null || criteria instanceof Error) {
    return 0
  }

  let result = 0;
  const isWildcard = criteria === '*';
  const tokenizedCriteria = isWildcard ? null : parse(criteria + '');

  for (let i = 0; i < range.length; i++) {
    const value = range[i];
    const sumValue = sum_range[i];

    if (isWildcard) {
      result += value;
    } else {
      const tokens = [createToken(value, TOKEN_TYPE_LITERAL)].concat(tokenizedCriteria);

      result += compute(tokens) ? sumValue : 0;
    }
  }

  return result
}

/**
 * Adds the values in a range that meet multiple criteria.
 *
 * Category: Math and trigonometry
 *
 * @returns
 */
function SUMIFS() {
  const args = argsToArray(arguments);
  const range = parseNumberArray(flatten(args.shift()));

  if (range instanceof Error) {
    return range
  }

  const criterias = args;
  const criteriaLength = criterias.length / 2;

  for (let i = 0; i < criteriaLength; i++) {
    criterias[i * 2] = flatten(criterias[i * 2]);
  }

  let result = 0;

  for (let i = 0; i < range.length; i++) {
    let isMeetCondition = false;

    for (let j = 0; j < criteriaLength; j++) {
      const valueToTest = criterias[j * 2][i];
      const criteria = criterias[j * 2 + 1];
      const isWildcard = criteria === void 0 || criteria === '*';
      let computedResult = false;

      if (isWildcard) {
        computedResult = true;
      } else {
        const tokenizedCriteria = parse(criteria + '');
        const tokens = [createToken(valueToTest, TOKEN_TYPE_LITERAL)].concat(
          tokenizedCriteria
        );

        computedResult = compute(tokens);
      }

      // Criterias are calculated as AND so any `false` breakes the loop as unmeet condition
      if (!computedResult) {
        isMeetCondition = false;
        break
      }

      isMeetCondition = true;
    }

    if (isMeetCondition) {
      result += range[i];
    }
  }

  return result
}

/**
 * Returns the sum of the products of corresponding array components.
 *
 * Category: Math and trigonometry
 *
 * @returns
 */
function SUMPRODUCT() {
  if (!arguments || arguments.length === 0) {
    return value
  }

  const arrays = arguments.length + 1;
  let result = 0;
  let product;
  let k;
  let _i;
  let _ij;

  for (let i = 0; i < arguments[0].length; i++) {
    if (!(arguments[0][i] instanceof Array)) {
      product = 1;

      for (k = 1; k < arrays; k++) {
        const _i_arg = arguments[k - 1][i];

        if (_i_arg instanceof Error) {
          return _i_arg
        }

        _i = parseNumber(_i_arg);

        if (_i instanceof Error) {
          return _i
        }

        product *= _i;
      }

      result += product;
    } else {
      for (let j = 0; j < arguments[0][i].length; j++) {
        product = 1;

        for (k = 1; k < arrays; k++) {
          const _ij_arg = arguments[k - 1][i][j];

          if (_ij_arg instanceof Error) {
            return _ij_arg
          }

          _ij = parseNumber(_ij_arg);

          if (_ij instanceof Error) {
            return _ij
          }

          product *= _ij;
        }

        result += product;
      }
    }
  }

  return result
}

/**
 * Returns the sum of the squares of the arguments.
 *
 * Category: Math and trigonometry
 *
 * @param {*} args number1, number2, ... Number1 is required, subsequent numbers are optional. 1 to 255 arguments for which you want the sum of the squares. You can also use a single array or a reference to an array instead of arguments separated by commas.
 * @returns
 */
function SUMSQ() {
  const numbers = parseNumberArray(flatten(arguments));

  if (numbers instanceof Error) {
    return numbers
  }

  let result = 0;
  const length = numbers.length;

  for (let i = 0; i < length; i++) {
    result += ISNUMBER(numbers[i]) ? numbers[i] * numbers[i] : 0;
  }

  return result
}

/**
 * Returns the sum of the difference of squares of corresponding values in two arrays.
 *
 * Category: Math and trigonometry
 *
 * @param {*} array_x The first array or range of values.
 * @param {*} array_y The second array or range of values.
 * @returns
 */
function SUMX2MY2(array_x, array_y) {
  array_x = parseNumberArray(flatten(array_x));
  array_y = parseNumberArray(flatten(array_y));

  if (anyIsError(array_x, array_y)) {
    return value
  }

  let result = 0;

  for (let i = 0; i < array_x.length; i++) {
    result += array_x[i] * array_x[i] - array_y[i] * array_y[i];
  }

  return result
}

/**
 * Returns the sum of the sum of squares of corresponding values in two arrays.
 *
 * Category: Math and trigonometry
 *
 * @param {*} array_x The first array or range of values.
 * @param {*} array_y The second array or range of values.
 * @returns
 */
function SUMX2PY2(array_x, array_y) {
  array_x = parseNumberArray(flatten(array_x));
  array_y = parseNumberArray(flatten(array_y));

  if (anyIsError(array_x, array_y)) {
    return value
  }

  let result = 0;
  array_x = parseNumberArray(flatten(array_x));
  array_y = parseNumberArray(flatten(array_y));

  for (let i = 0; i < array_x.length; i++) {
    result += array_x[i] * array_x[i] + array_y[i] * array_y[i];
  }

  return result
}

/**
 * Returns the sum of squares of differences of corresponding values in two arrays.
 *
 * Category: Math and trigonometry
 *
 * @param {*} array_x The first array or range of values.
 * @param {*} array_y The second array or range of values.
 * @returns
 */
function SUMXMY2(array_x, array_y) {
  array_x = parseNumberArray(flatten(array_x));
  array_y = parseNumberArray(flatten(array_y));

  if (anyIsError(array_x, array_y)) {
    return value
  }

  let result = 0;
  array_x = flatten(array_x);
  array_y = flatten(array_y);

  for (let i = 0; i < array_x.length; i++) {
    result += Math.pow(array_x[i] - array_y[i], 2);
  }

  return result
}

/**
 * Returns the tangent of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The angle in radians for which you want the tangent.
 * @returns
 */
function TAN(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  return Math.tan(number)
}

/**
 * Returns the hyperbolic tangent of a number.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number Any real number.
 * @returns
 */
function TANH(number) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  const e2 = Math.exp(2 * number);

  return (e2 - 1) / (e2 + 1)
}

/**
 * Truncates a number to an integer.
 *
 * Category: Math and trigonometry
 *
 * @param {*} number The number you want to truncate.
 * @param {*} num_digits Optional. A number specifying the precision of the truncation. The default value for num_digits is 0 (zero).
 * @returns
 */
function TRUNC(number, num_digits) {
  number = parseNumber(number);
  num_digits = parseNumber(num_digits);
  const anyError$1 = anyError(number, num_digits);

  if (anyError$1) {
    return anyError$1
  }

  const sign = number > 0 ? 1 : -1;

  return (sign * Math.floor(Math.abs(number) * Math.pow(10, num_digits))) / Math.pow(10, num_digits)
}

/**
 * Following functions are part of Formula.js only and not found in Excel.
 * Purpose of these functions is to implement an Excel like behaviour for fundamental math symbols such as +, -, /, *, etc.
 */

/**
 * Formula.js only.
 *
 * @param {*} num1
 * @param {*} num2
 * @returns
 */
function ADD(num1, num2) {
  if (arguments.length !== 2) {
    return na
  }

  num1 = parseNumber(num1);
  num2 = parseNumber(num2);
  const anyError$1 = anyError(num1, num2);

  if (anyError$1) {
    return anyError$1
  }

  return num1 + num2
}

/**
 * Formula.js only
 *
 * @param {*} dividend
 * @param {*} divisor
 * @returns
 */
function DIVIDE(dividend, divisor) {
  if (arguments.length !== 2) {
    return na
  }

  dividend = parseNumber(dividend);
  divisor = parseNumber(divisor);
  const anyError$1 = anyError(dividend, divisor);

  if (anyError$1) {
    return anyError$1
  }

  if (divisor === 0) {
    return div0
  }

  return dividend / divisor
}

/**
 * Formula.js only
 *
 * @param {*} value1
 * @param {*} value2
 * @returns
 */
function EQ(value1, value2) {
  if (arguments.length !== 2) {
    return na
  }

  if (value1 instanceof Error) {
    return value1
  }

  if (value2 instanceof Error) {
    return value2
  }

  if (value1 === null) {
    value1 = undefined;
  }

  if (value2 === null) {
    value2 = undefined;
  }

  return value1 === value2
}

/**
 * Formula.js only
 *
 * @param {*} num1
 * @param {*} num2
 * @returns
 */
function GT(num1, num2) {
  if (arguments.length !== 2) {
    return na
  }

  if (num1 instanceof Error) {
    return num1
  }

  if (num2 instanceof Error) {
    return num2
  }

  if (anyIsString(num1, num2)) {
    num1 = parseString(num1);
    num2 = parseString(num2);
  } else {
    num1 = parseNumber(num1);
    num2 = parseNumber(num2);
  }

  const anyError$1 = anyError(num1, num2);

  if (anyError$1) {
    return anyError$1
  }

  return num1 > num2
}

/**
 * Formula.js only
 *
 * @param {*} num1
 * @param {*} num2
 * @returns
 */
function GTE(num1, num2) {
  if (arguments.length !== 2) {
    return na
  }

  if (anyIsString(num1, num2)) {
    num1 = parseString(num1);
    num2 = parseString(num2);
  } else {
    num1 = parseNumber(num1);
    num2 = parseNumber(num2);
  }

  const anyError$1 = anyError(num1, num2);

  if (anyError$1) {
    return anyError$1
  }

  return num1 >= num2
}

/**
 * Formula.js only
 *
 * @param {*} num1
 * @param {*} num2
 * @returns
 */
function LT(num1, num2) {
  if (arguments.length !== 2) {
    return na
  }

  if (anyIsString(num1, num2)) {
    num1 = parseString(num1);
    num2 = parseString(num2);
  } else {
    num1 = parseNumber(num1);
    num2 = parseNumber(num2);
  }

  const anyError$1 = anyError(num1, num2);

  if (anyError$1) {
    return anyError$1
  }

  return num1 < num2
}

/**
 * Formula.js only
 *
 * @param {*} num1
 * @param {*} num2
 * @returns
 */
function LTE(num1, num2) {
  if (arguments.length !== 2) {
    return na
  }

  if (anyIsString(num1, num2)) {
    num1 = parseString(num1);
    num2 = parseString(num2);
  } else {
    num1 = parseNumber(num1);
    num2 = parseNumber(num2);
  }

  const anyError$1 = anyError(num1, num2);

  if (anyError$1) {
    return anyError$1
  }

  return num1 <= num2
}

/**
 * Formula.js only
 *
 * @param {*} num1
 * @param {*} num2
 * @returns
 */
function MINUS(num1, num2) {
  if (arguments.length !== 2) {
    return na
  }

  num1 = parseNumber(num1);
  num2 = parseNumber(num2);
  const anyError$1 = anyError(num1, num2);

  if (anyError$1) {
    return anyError$1
  }

  return num1 - num2
}

/**
 * Formula.js only
 *
 * @param {*} factor1
 * @param {*} factor2
 * @returns
 */
function MULTIPLY(factor1, factor2) {
  if (arguments.length !== 2) {
    return na
  }

  factor1 = parseNumber(factor1);
  factor2 = parseNumber(factor2);
  const anyError$1 = anyError(factor1, factor2);

  if (anyError$1) {
    return anyError$1
  }

  return factor1 * factor2
}

/**
 * Formula.js only
 *
 * @param {*} value1
 * @param {*} value2
 * @returns
 */
function NE(value1, value2) {
  if (arguments.length !== 2) {
    return na
  }

  if (value1 instanceof Error) {
    return value1
  }

  if (value2 instanceof Error) {
    return value2
  }

  if (value1 === null) {
    value1 = undefined;
  }

  if (value2 === null) {
    value2 = undefined;
  }

  return value1 !== value2
}

/**
 * Formula.js only
 *
 * @param {*} base
 * @param {*} exponent
 * @returns
 */
function POW(base, exponent) {
  if (arguments.length !== 2) {
    return na
  }

  return POWER(base, exponent)
}

var symbols = /*#__PURE__*/Object.freeze({
  __proto__: null,
  ADD: ADD,
  DIVIDE: DIVIDE,
  EQ: EQ,
  GT: GT,
  GTE: GTE,
  LT: LT,
  LTE: LTE,
  MINUS: MINUS,
  MULTIPLY: MULTIPLY,
  NE: NE,
  POW: POW
});

const d1900 = new Date(Date.UTC(1900, 0, 1));
const WEEK_STARTS = [
  undefined,
  0,
  1,
  undefined,
  undefined,
  undefined,
  undefined,
  undefined,
  undefined,
  undefined,
  undefined,
  undefined,
  1,
  2,
  3,
  4,
  5,
  6,
  0
];
const WEEK_TYPES = [
  [],
  [1, 2, 3, 4, 5, 6, 7],
  [7, 1, 2, 3, 4, 5, 6],
  [6, 0, 1, 2, 3, 4, 5],
  [],
  [],
  [],
  [],
  [],
  [],
  [],
  [7, 1, 2, 3, 4, 5, 6],
  [6, 7, 1, 2, 3, 4, 5],
  [5, 6, 7, 1, 2, 3, 4],
  [4, 5, 6, 7, 1, 2, 3],
  [3, 4, 5, 6, 7, 1, 2],
  [2, 3, 4, 5, 6, 7, 1],
  [1, 2, 3, 4, 5, 6, 7]
];
const WEEKEND_TYPES = [
  [],
  [6, 0],
  [0, 1],
  [1, 2],
  [2, 3],
  [3, 4],
  [4, 5],
  [5, 6],
  undefined,
  undefined,
  undefined,
  [0, 0],
  [1, 1],
  [2, 2],
  [3, 3],
  [4, 4],
  [5, 5],
  [6, 6]
];

/**
 * Returns the serial number of a particular date.
 *
 * Category: Date and time
 *
 * @param {*} year Year
 * @param {*} month Month
 * @param {*} day Day
 * @returns
 */
function DATE(year, month, day) {
  let result;

  year = parseNumber(year);
  month = parseNumber(month);
  day = parseNumber(day);

  if (anyIsError(year, month, day)) {
    result = value;
  } else {
    result = new Date(year, month - 1, day);

    if (result.getFullYear() < 0) {
      result = num;
    }
  }

  return result
}

/**
 * Calculates the number of days, months, or years between two dates. This function is useful in formulas where you need to calculate an age.
 *
 * Category: Date and time
 *
 * @param {*} start_date A date that represents the first, or starting date of a given period.
 * @param {*} end_date A date that represents the last, or ending, date of the period.
 * @param {*} unit The type of information that you want returned, where:
 - "Y": The number of complete years in the period.
 - "M": The number of complete months in the period.
 - "D": The number of days in the period.
 - "MD": The difference between the days in start_date and end_date. The months and years of the dates are ignored.
 - "YM": The difference between the months in start_date and end_date. The days and years of the dates are ignored
 - "YD": The difference between the days of start_date and end_date. The years of the dates are ignored.
 * @returns
 */
function DATEDIF(start_date, end_date, unit) {
  unit = unit.toUpperCase();
  start_date = parseDate(start_date);
  end_date = parseDate(end_date);

  const start_date_year = start_date.getFullYear();
  const start_date_month = start_date.getMonth();
  const start_date_day = start_date.getDate();
  const end_date_year = end_date.getFullYear();
  const end_date_month = end_date.getMonth();
  const end_date_day = end_date.getDate();

  let result;

  switch (unit) {
    case 'Y':
      result = Math.floor(YEARFRAC(start_date, end_date));
      break
    case 'D':
      result = DAYS(end_date, start_date);
      break
    case 'M':
      result = end_date_month - start_date_month + 12 * (end_date_year - start_date_year);

      if (end_date_day < start_date_day) {
        result--;
      }

      break
    case 'MD':
      if (start_date_day <= end_date_day) {
        result = end_date_day - start_date_day;
      } else {
        if (end_date_month === 0) {
          start_date.setFullYear(end_date_year - 1);
          start_date.setMonth(12);
        } else {
          start_date.setFullYear(end_date_year);
          start_date.setMonth(end_date_month - 1);
        }

        result = DAYS(end_date, start_date);
      }

      break
    case 'YM':
      result = end_date_month - start_date_month + 12 * (end_date_year - start_date_year);

      if (end_date_day < start_date_day) {
        result--;
      }

      result = result % 12;
      break
    case 'YD':
      if (end_date_month > start_date_month || (end_date_month === start_date_month && end_date_day < start_date_day)) {
        start_date.setFullYear(end_date_year);
      } else {
        start_date.setFullYear(end_date_year - 1);
      }

      result = DAYS(end_date, start_date);
      break
  }

  return result
}

/**
 * Converts a date in the form of text to a serial number.
 *
 * Category: Date and time
 *
 * @param {*} date_text Text that represents a date in an Excel date format, or a reference to a value that contains text that represents a date in an Excel date format.
 * @returns
 */
function DATEVALUE(date_text) {
  if (typeof date_text !== 'string') {
    return value
  }

  const date = Date.parse(date_text);

  if (isNaN(date)) {
    return value
  }

  return new Date(date_text)
}

/**
 * Converts a serial number to a day of the month.
 *
 * Category: Date and time
 *
 * @param {*} serial_number The date of the day you are trying to find.
 * @returns
 */
function DAY(serial_number) {
  const date = parseDate(serial_number);

  if (date instanceof Error) {
    return date
  }

  return date.getDate()
}

function startOfDay(date) {
  const newDate = new Date(date);
  newDate.setHours(0, 0, 0, 0);

  return newDate
}

/**
 * Returns the number of days between two dates.
 *
 * Category: Date and time
 *
 * @param {*} end_date Start_date and End_date are the two dates between which you want to know the number of days.
 * @param {*} start_date Start_date and End_date are the two dates between which you want to know the number of days.
 * @returns
 */
function DAYS(end_date, start_date) {
  end_date = parseDate(end_date);
  start_date = parseDate(start_date);

  if (end_date instanceof Error) {
    return end_date
  }

  if (start_date instanceof Error) {
    return start_date
  }

  return serial(startOfDay(end_date)) - serial(startOfDay(start_date))
}

/**
 * Calculates the number of days between two dates based on a 360-day year.
 *
 * Category: Date and time
 *
 * @param {*} start_date A date that represents the start date. If start_date occurs after end_date, the DAYS360 function returns a negative number.
 * @param {*} end_date A date that represents the end date.
 * @param {*} method Optional. A logical value that specifies whether to use the U.S. or European method in the calculation.
 * @returns
 */
function DAYS360(start_date, end_date, method) {
  method = parseBool(method || 'false');
  start_date = parseDate(start_date);
  end_date = parseDate(end_date);

  if (start_date instanceof Error) {
    return start_date
  }

  if (end_date instanceof Error) {
    return end_date
  }

  if (method instanceof Error) {
    return method
  }

  const sm = start_date.getMonth();
  let em = end_date.getMonth();
  let sd, ed;

  if (method) {
    sd = start_date.getDate() === 31 ? 30 : start_date.getDate();
    ed = end_date.getDate() === 31 ? 30 : end_date.getDate();
  } else {
    const smd = new Date(start_date.getFullYear(), sm + 1, 0).getDate();
    const emd = new Date(end_date.getFullYear(), em + 1, 0).getDate();
    sd = start_date.getDate() === smd ? 30 : start_date.getDate();

    if (end_date.getDate() === emd) {
      if (sd < 30) {
        em++;
        ed = 1;
      } else {
        ed = 30;
      }
    } else {
      ed = end_date.getDate();
    }
  }

  return 360 * (end_date.getFullYear() - start_date.getFullYear()) + 30 * (em - sm) + (ed - sd)
}

/**
 * Returns the serial number of the date that is the indicated number of months before or after the start date.
 *
 * Category: Date and time
 *
 * @param {*} start_date A date that represents the start date.
 * @param {*} months The number of months before or after start_date. A positive value for months yields a future date; a negative value yields a past date.
 * @returns
 */
function EDATE(start_date, months) {
  start_date = parseDate(start_date);

  if (start_date instanceof Error) {
    return start_date
  }

  if (isNaN(months)) {
    return value
  }

  months = parseInt(months, 10);
  start_date.setMonth(start_date.getMonth() + months);

  return start_date
}

/**
 * Returns the serial number of the last day of the month before or after a specified number of months.
 *
 * Category: Date and time
 *
 * @param {*} start_date A date that represents the starting date.
 * @param {*} months The number of months before or after start_date. A positive value for months yields a future date; a negative value yields a past date.
 * @returns
 */
function EOMONTH(start_date, months) {
  start_date = parseDate(start_date);

  if (start_date instanceof Error) {
    return start_date
  }

  if (isNaN(months)) {
    return value
  }

  months = parseInt(months, 10);

  return new Date(start_date.getFullYear(), start_date.getMonth() + months + 1, 0)
}

/**
 * Converts a serial number to an hour.
 *
 * Category: Date and time
 *
 * @param {*} serial_number The time that contains the hour you want to find. Times may be entered as text strings within quotation marks (for example, "6:45 PM"), as decimal numbers (for example, 0.78125, which represents 6:45 PM), or as results of other formulas or functions (for example, TIMEVALUE("6:45 PM")).
 * @returns
 */
function HOUR(serial_number) {
  serial_number = parseDate(serial_number);

  if (serial_number instanceof Error) {
    return serial_number
  }

  return serial_number.getHours()
}

/**
 * Returns the number of the ISO week number of the year for a given date.
 *
 * Category: Date and time
 *
 * @param {*} date Date is the date-time code used by Excel for date and time calculation.
 * @returns
 */
function ISOWEEKNUM(date) {
  date = parseDate(date);

  if (date instanceof Error) {
    return date
  }

  date = startOfDay(date);
  date.setDate(date.getDate() + 4 - (date.getDay() || 7));
  const yearStart = new Date(date.getFullYear(), 0, 1);

  return Math.ceil(((date - yearStart) / 86400000 + 1) / 7)
}

/**
 * Converts a serial number to a minute.
 *
 * Category: Date and time
 *
 * @param {*} serial_number The time that contains the minute you want to find. Times may be entered as text strings within quotation marks (for example, "6:45 PM"), as decimal numbers (for example, 0.78125, which represents 6:45 PM), or as results of other formulas or functions (for example, TIMEVALUE("6:45 PM")).
 * @returns
 */
function MINUTE(serial_number) {
  serial_number = parseDate(serial_number);

  if (serial_number instanceof Error) {
    return serial_number
  }

  return serial_number.getMinutes()
}

/**
 * Converts a serial number to a month.
 *
 * Category: Date and time
 *
 * @param {*} serial_number The date of the month you are trying to find.
 * @returns
 */
function MONTH(serial_number) {
  serial_number = parseDate(serial_number);

  if (serial_number instanceof Error) {
    return serial_number
  }

  return serial_number.getMonth() + 1
}

/**
 * Returns the number of whole workdays between two dates.
 *
 * Category: Date and time
 *
 * @param {*} start_date A date that represents the start date.
 * @param {*} end_date A date that represents the end date.
 * @param {*} holidays Optional. An optional range of one or more dates to exclude from the working calendar, such as state and federal holidays and floating holidays. The list can be either a range of values that contains the dates or an array constant of the serial numbers that represent the dates.
 * @returns
 */
function NETWORKDAYS(start_date, end_date, holidays) {
  return NETWORKDAYS.INTL(start_date, end_date, 1, holidays)
}

/**
 * Returns the number of whole workdays between two dates using parameters to indicate which and how many days are weekend days.
 *
 * Category: Date and time
 *
 * @param {*} start_date The date for from which the difference is to be computed. The start_date can be earlier than, the same as, or later than the end_date.
 * @param {*} end_date The date for to which the difference is to be computed.
 * @param {*} weekend Optional. Indicates the days of the week that are weekend days and are not included in the number of whole working days between start_date and end_date. Weekend is a weekend number or string that specifies when weekends occur. Weekend number values indicate the following weekend days:
 * @param {*} holidays Optional. An optional set of one or more dates that are to be excluded from the working day calendar. holidays shall be a range of values that contain the dates, or an array constant of the serial values that represent those dates. The ordering of dates or serial values in holidays can be arbitrary.
 * @returns
 */
NETWORKDAYS.INTL = (start_date, end_date, weekend, holidays) => {
  start_date = parseDate(start_date);

  if (start_date instanceof Error) {
    return start_date
  }

  end_date = parseDate(end_date);

  if (end_date instanceof Error) {
    return end_date
  }

  let isMask = false;
  const maskDays = [];
  const maskIndex = [1, 2, 3, 4, 5, 6, 0];
  const maskRegex = new RegExp('^[0|1]{7}$');

  if (weekend === undefined) {
    weekend = WEEKEND_TYPES[1];
  } else if (typeof weekend === 'string' && maskRegex.test(weekend)) {
    isMask = true;
    weekend = weekend.split('');

    for (let i = 0; i < weekend.length; i++) {
      if (weekend[i] === '1') {
        maskDays.push(maskIndex[i]);
      }
    }
  } else {
    weekend = WEEKEND_TYPES[weekend];
  }

  if (!(weekend instanceof Array)) {
    return value
  }

  if (holidays === undefined) {
    holidays = [];
  } else if (!(holidays instanceof Array)) {
    holidays = [holidays];
  }

  for (let i = 0; i < holidays.length; i++) {
    const h = parseDate(holidays[i]);

    if (h instanceof Error) {
      return h
    }

    holidays[i] = h;
  }

  const days = Math.round((end_date - start_date) / (1000 * 60 * 60 * 24)) + 1;
  let total = days;
  const day = start_date;

  for (let i = 0; i < days; i++) {
    const d = new Date().getTimezoneOffset() > 0 ? day.getUTCDay() : day.getDay();
    let dec = isMask ? maskDays.includes(d) : d === weekend[0] || d === weekend[1];

    for (let j = 0; j < holidays.length; j++) {
      const holiday = holidays[j];

      if (
        holiday.getDate() === day.getDate() &&
        holiday.getMonth() === day.getMonth() &&
        holiday.getFullYear() === day.getFullYear()
      ) {
        dec = true;
        break
      }
    }

    if (dec) {
      total--;
    }

    day.setDate(day.getDate() + 1);
  }

  return total
};

/**
 * Returns the serial number of the current date and time.
 *
 * Category: Date and time
 *
 * @returns
 */
function NOW() {
  return new Date()
}

/**
 * Converts a serial number to a second.
 *
 * Category: Date and time
 *
 * @param {*} serial_number The time that contains the seconds you want to find. Times may be entered as text strings within quotation marks (for example, "6:45 PM"), as decimal numbers (for example, 0.78125, which represents 6:45 PM), or as results of other formulas or functions (for example, TIMEVALUE("6:45 PM")).
 * @returns
 */
function SECOND(serial_number) {
  serial_number = parseDate(serial_number);

  if (serial_number instanceof Error) {
    return serial_number
  }

  return serial_number.getSeconds()
}

/**
 * Returns the serial number of a particular time.
 *
 * Category: Date and time
 *
 * @param {*} hour A number from 0 (zero) to 32767 representing the hour. Any value greater than 23 will be divided by 24 and the remainder will be treated as the hour value. For example, TIME(27,0,0) = TIME(3,0,0) = .125 or 3:00 AM.
 * @param {*} minute A number from 0 to 32767 representing the minute. Any value greater than 59 will be converted to hours and minutes. For example, TIME(0,750,0) = TIME(12,30,0) = .520833 or 12:30 PM.
 * @param {*} second A number from 0 to 32767 representing the second. Any value greater than 59 will be converted to hours, minutes, and seconds. For example, TIME(0,0,2000) = TIME(0,33,22) = .023148 or 12:33:20 AM
 * @returns
 */
function TIME(hour, minute, second) {
  hour = parseNumber(hour);
  minute = parseNumber(minute);
  second = parseNumber(second);

  if (anyIsError(hour, minute, second)) {
    return value
  }

  if (hour < 0 || minute < 0 || second < 0) {
    return num
  }

  return (3600 * hour + 60 * minute + second) / 86400
}

/**
 * Converts a time in the form of text to a serial number.
 *
 * Category: Date and time
 *
 * @param {*} time_text A text string that represents a time in any one of the Microsoft Excel time formats; for example, "6:45 PM" and "18:45" text strings within quotation marks that represent time.
 * @returns
 */
function TIMEVALUE(time_text) {
  time_text = parseDate(time_text);

  if (time_text instanceof Error) {
    return time_text
  }

  return (3600 * time_text.getHours() + 60 * time_text.getMinutes() + time_text.getSeconds()) / 86400
}

/**
 * Returns the serial number of today's date.
 *
 * Category: Date and time
 *
 * @returns
 */
function TODAY() {
  return startOfDay(new Date())
}

/**
 * Converts a serial number to a day of the week.
 *
 * Category: Date and time
 *
 * @param {*} serial_number A sequential number that represents the date of the day you are trying to find.
 * @param {*} return_type Optional. A number that determines the type of return value.
 * @returns
 */
function WEEKDAY(serial_number, return_type) {
  serial_number = parseDate(serial_number);

  if (serial_number instanceof Error) {
    return serial_number
  }

  if (return_type === undefined) {
    return_type = 1;
  }

  const day = serial_number.getDay();

  return WEEK_TYPES[return_type][day]
}

/**
 * Converts a serial number to a number representing where the week falls numerically with a year.
 *
 * Category: Date and time
 *
 * @param {*} serial_number A date within the week.
 * @param {*} return_type Optional. A number that determines on which day the week begins. The default is 1.
 * @returns
 */
function WEEKNUM(serial_number, return_type) {
  serial_number = parseDate(serial_number);

  if (serial_number instanceof Error) {
    return serial_number
  }

  if (return_type === undefined) {
    return_type = 1;
  }

  if (return_type === 21) {
    return ISOWEEKNUM(serial_number)
  }

  const week_start = WEEK_STARTS[return_type];
  let jan = new Date(serial_number.getFullYear(), 0, 1);
  const inc = jan.getDay() < week_start ? 1 : 0;
  jan -= Math.abs(jan.getDay() - week_start) * 24 * 60 * 60 * 1000;

  return Math.floor((serial_number - jan) / (1000 * 60 * 60 * 24) / 7 + 1) + inc
}

/**
 * Returns the serial number of the date before or after a specified number of workdays.
 *
 * Category: Date and time
 *
 * @param {*} start_date A date that represents the start date.
 * @param {*} days The number of nonweekend and nonholiday days before or after start_date. A positive value for days yields a future date; a negative value yields a past date.
 * @param {*} holidays Optional. An optional list of one or more dates to exclude from the working calendar, such as state and federal holidays and floating holidays. The list can be either a range of values that contain the dates or an array constant of the serial numbers that represent the dates.
 * @returns
 */
function WORKDAY(start_date, days, holidays) {
  return WORKDAY.INTL(start_date, days, 1, holidays)
}

/**
 * Returns the serial number of the date before or after a specified number of workdays using parameters to indicate which and how many days are weekend days.
 *
 * Category: Date and time
 *
 * @param {*} start_date The start date, truncated to integer.
 * @param {*} days The number of workdays before or after the start_date. A positive value yields a future date; a negative value yields a past date; a zero value yields the start_date. Day-offset is truncated to an integer.
 * @param {*} weekend Optional. Indicates the days of the week that are weekend days and are not considered working days. Weekend is a weekend number or string that specifies when weekends occur. Weekend number values indicate the following weekend days:
 * @param {*} holidays Optional. An optional set of one or more dates that are to be excluded from the working day calendar. Holidays shall be a range of values that contain the dates, or an array constant of the serial values that represent those dates. The ordering of dates or serial values in holidays can be arbitrary.
 * @returns
 */
WORKDAY.INTL = (start_date, days, weekend, holidays) => {
  start_date = parseDate(start_date);

  if (start_date instanceof Error) {
    return start_date
  }

  days = parseNumber(days);

  if (days instanceof Error) {
    return days
  }

  if (days < 0) {
    return num
  }

  if (weekend === undefined) {
    weekend = WEEKEND_TYPES[1];
  } else {
    weekend = WEEKEND_TYPES[weekend];
  }

  if (!(weekend instanceof Array)) {
    return value
  }

  if (holidays === undefined) {
    holidays = [];
  } else if (!(holidays instanceof Array)) {
    holidays = [holidays];
  }

  for (let i = 0; i < holidays.length; i++) {
    const h = parseDate(holidays[i]);

    if (h instanceof Error) {
      return h
    }

    holidays[i] = h;
  }

  let d = 0;

  while (d < days) {
    start_date.setDate(start_date.getDate() + 1);
    const day = start_date.getDay();

    if (day === weekend[0] || day === weekend[1]) {
      continue
    }

    for (let j = 0; j < holidays.length; j++) {
      const holiday = holidays[j];

      if (
        holiday.getDate() === start_date.getDate() &&
        holiday.getMonth() === start_date.getMonth() &&
        holiday.getFullYear() === start_date.getFullYear()
      ) {
        d--;
        break
      }
    }

    d++;
  }

  return start_date
};

/**
 * Converts a serial number to a year.
 *
 * Category: Date and time
 *
 * @param {*} serial_number The date of the year you want to find.
 * @returns
 */
function YEAR(serial_number) {
  serial_number = parseDate(serial_number);

  if (serial_number instanceof Error) {
    return serial_number
  }

  return serial_number.getFullYear()
}

function isLeapYear(year) {
  return new Date(year, 1, 29).getMonth() === 1
}

// TODO : Use DAYS ?
function daysBetween(start_date, end_date) {
  return Math.ceil((end_date - start_date) / 1000 / 60 / 60 / 24)
}

/**
 * Returns the year fraction representing the number of whole days between start_date and end_date.
 *
 * Category: Date and time
 *
 * @param {*} start_date A date that represents the start date.
 * @param {*} end_date A date that represents the end date.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function YEARFRAC(start_date, end_date, basis) {
  start_date = parseDate(start_date);

  if (start_date instanceof Error) {
    return start_date
  }

  end_date = parseDate(end_date);

  if (end_date instanceof Error) {
    return end_date
  }

  basis = basis || 0;
  let sd = start_date.getDate();
  const sm = start_date.getMonth() + 1;
  const sy = start_date.getFullYear();
  let ed = end_date.getDate();
  const em = end_date.getMonth() + 1;
  const ey = end_date.getFullYear();

  switch (basis) {
    case 0:
      // US (NASD) 30/360
      if (sd === 31 && ed === 31) {
        sd = 30;
        ed = 30;
      } else if (sd === 31) {
        sd = 30;
      } else if (sd === 30 && ed === 31) {
        ed = 30;
      }

      return (ed + em * 30 + ey * 360 - (sd + sm * 30 + sy * 360)) / 360
    case 1: {
      // Actual/actual
      const feb29Between = (date1, date2) => {
        const year1 = date1.getFullYear();
        const mar1year1 = new Date(year1, 2, 1);

        if (isLeapYear(year1) && date1 < mar1year1 && date2 >= mar1year1) {
          return true
        }

        const year2 = date2.getFullYear();
        const mar1year2 = new Date(year2, 2, 1);

        return isLeapYear(year2) && date2 >= mar1year2 && date1 < mar1year2
      };

      let ylength = 365;

      if (sy === ey || (sy + 1 === ey && (sm > em || (sm === em && sd >= ed)))) {
        if ((sy === ey && isLeapYear(sy)) || feb29Between(start_date, end_date) || (em === 1 && ed === 29)) {
          ylength = 366;
        }

        return daysBetween(start_date, end_date) / ylength
      }

      const years = ey - sy + 1;
      const days = (new Date(ey + 1, 0, 1) - new Date(sy, 0, 1)) / 1000 / 60 / 60 / 24;
      const average = days / years;

      return daysBetween(start_date, end_date) / average
    }

    case 2:
      // Actual/360

      return daysBetween(start_date, end_date) / 360
    case 3:
      // Actual/365

      return daysBetween(start_date, end_date) / 365
    case 4:
      // European 30/360

      return (ed + em * 30 + ey * 360 - (sd + sm * 30 + sy * 360)) / 360
  }
}

function serial(date) {
  const addOn = date > -2203891200000 ? 2 : 1;

  return Math.ceil((date - d1900) / 86400000) + addOn
}

// TODO
/**
 * -- Not implemented --
 *
 * Changes full-width (double-byte) English letters or katakana within a character string to half-width (single-byte) characters.
 *
 * Category: Text
 *
 * @param {*} text The text or a reference to a value that contains the text you want to change. If text does not contain any full-width letters, text is not changed.
 * @returns
 */
function ASC() {
  throw new Error('ASC is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Converts a number to text, using the ß (baht) currency format.
 *
 * Category: Text
 *
 * @param {*} number A number you want to convert to text, or a reference to a value containing a number, or a formula that evaluates to a number.
 * @returns
 */
function BAHTTEXT() {
  throw new Error('BAHTTEXT is not implemented')
}

/**
 * Returns the character specified by the code number.
 *
 * Category: Text
 *
 * @param {*} number A number between 1 and 255 specifying which character you want. The character is from the character set used by your computer. Note: Excel for the web supports only CHAR(9), CHAR(10), CHAR(13), and CHAR(32) and above.
 * @returns
 */
function CHAR(number) {
  number = parseNumber(number);

  if (number === 0) {
    return value
  }

  if (number instanceof Error) {
    return number
  }

  return String.fromCharCode(number)
}

/**
 * Removes all nonprintable characters from text.
 *
 * Category: Text
 *
 * @param {*} text Any worksheet information from which you want to remove nonprintable characters.
 * @returns
 */
function CLEAN(text) {
  if (anyIsError(text)) {
    return text
  }

  text = text || '';
  const re = /[\0-\x1F]/g;

  return text.replace(re, '')
}

/**
 * Returns a numeric code for the first character in a text string.
 *
 * Category: Text
 *
 * @param {*} text The text for which you want the code of the first character.
 * @returns
 */
function CODE(text) {
  if (anyIsError(text)) {
    return text
  }

  text = text || '';
  let result = text.charCodeAt(0);

  if (isNaN(result)) {
    result = value;
  }

  return result
}

/**
 * Joins several text items into one text item.
 *
 * Category: Text
 *
 * @returns
 */
function CONCATENATE() {
  const args = flatten(arguments);
  const someError = anyError.apply(undefined, args);

  if (someError) {
    return someError
  }

  let trueFound = 0;

  while ((trueFound = args.indexOf(true)) > -1) {
    args[trueFound] = 'TRUE';
  }

  let falseFound = 0;

  while ((falseFound = args.indexOf(false)) > -1) {
    args[falseFound] = 'FALSE';
  }

  return args.join('')
}

const CONCAT = CONCATENATE;

// TODO
/**
 * -- Not implemented --
 *
 * Changes half-width (single-byte) English letters or katakana within a character string to full-width (double-byte) characters.
 *
 * Category: Text
 *
 * @param {*} text The text or a reference to a value that contains the text you want to change. If text does not contain any half-width English letters or katakana, text is not changed.
 * @returns
 */
function DBCS() {
  throw new Error('DBCS is not implemented')
}

/**
 * Converts a number to text, using the $ (dollar) currency format.
 *
 * Category: Text
 *
 * @param {*} number A number, a reference to a value containing a number, or a formula that evaluates to a number.
 * @param {*} decimals Optional. The number of digits to the right of the decimal point. If this is negative, the number is rounded to the left of the decimal point. If you omit decimals, it is assumed to be 2.
 * @returns
 */
function DOLLAR(number, decimals = 2) {
  number = parseNumber(number);
  if (isNaN(number)) {
    return value
  }

  number = ROUND(number, decimals);

  const options = {
    style: 'currency',
    currency: 'USD',
    minimumFractionDigits: decimals >= 0 ? decimals : 0,
    maximumFractionDigits: decimals >= 0 ? decimals : 0
  };

  const formattedNumber = number.toLocaleString('en-US', options);

  if (number < 0) {
    return '$(' + formattedNumber.slice(2) + ')'
  }

  return formattedNumber
}

/**
 * Checks to see if two text values are identical.
 *
 * Category: Text
 *
 * @param {*} text1 The first text string.
 * @param {*} text2 The second text string.
 * @returns
 */
function EXACT(text1, text2) {
  if (arguments.length !== 2) {
    return na
  }

  const someError = anyError(text1, text2);

  if (someError) {
    return someError
  }

  text1 = parseString(text1);
  text2 = parseString(text2);

  return text1 === text2
}

/**
 * Locate one text string within a second text string, and return the number of the starting position of the first text string from the first character of the second text string.
 *
 * Category: Text
 *
 * @param {*} find_text The text you want to find.
 * @param {*} within_text The text containing the text you want to find.
 * @param {*} start_num Optional. Specifies the character at which to start the search. The first character in within_text is character number 1. If you omit start_num, it is assumed to be 1.
 * @returns
 */
function FIND(find_text, within_text, start_num) {
  if (arguments.length < 2) {
    return na
  }

  find_text = parseString(find_text);
  within_text = parseString(within_text);
  start_num = start_num === undefined ? 0 : start_num;
  const found_index = within_text.indexOf(find_text, start_num - 1);

  if (found_index === -1) {
    return value
  }

  return found_index + 1
}

/**
 * Formats a number as text with a fixed number of decimals.
 *
 * Category: Text
 *
 * @param {*} number The number you want to round and convert to text.
 * @param {*} decimals Optional. The number of digits to the right of the decimal point.
 * @param {*} no_commas Optional. A logical value that, if TRUE, prevents FIXED from including commas in the returned text.
 * @returns
 */
function FIXED(number, decimals = 2, no_commas = false) {
  number = parseNumber(number);
  if (isNaN(number)) {
    return value
  }

  decimals = parseNumber(decimals);
  if (isNaN(decimals)) {
    return value
  }

  if (decimals < 0) {
    const factor = Math.pow(10, -decimals);
    number = Math.round(number / factor) * factor;
  } else {
    number = number.toFixed(decimals);
  }

  if (no_commas) {
    number = number.toString().replace(/,/g, '');
  } else {
    const parts = number.toString().split('.');
    parts[0] = parts[0].replace(/\B(?=(\d{3})+$)/g, ',');
    number = parts.join('.');
  }

  return number
}

/**
 * Returns the leftmost characters from a text value.
 *
 * Category: Text
 *
 * @param {*} text The text string that contains the characters you want to extract.
 * @param {*} num_chars Optional. Specifies the number of characters you want LEFT to extract.
 * @returns
 */
function LEFT(text, num_chars) {
  const someError = anyError(text, num_chars);

  if (someError) {
    return someError
  }

  text = parseString(text);
  num_chars = num_chars === undefined ? 1 : num_chars;
  num_chars = parseNumber(num_chars);

  if (num_chars instanceof Error || typeof text !== 'string') {
    return value
  }

  return text.substring(0, num_chars)
}

/**
 * Returns the number of characters in a text string
 *
 * Category: Text
 *
 * @param {*} text The text whose length you want to find. Spaces count as characters.
 * @returns
 */
function LEN(text) {
  if (arguments.length === 0) {
    return error
  }

  if (text instanceof Error) {
    return text
  }

  if (Array.isArray(text)) {
    return value
  }

  const textAsString = parseString(text);

  return textAsString.length
}

/**
 * Converts text to lowercase.
 *
 * Category: Text
 *
 * @param {*} text The text you want to convert to lowercase. LOWER does not change characters in text that are not letters.
 * @returns
 */
function LOWER(text) {
  if (arguments.length !== 1) {
    return value
  }

  text = parseString(text);

  if (anyIsError(text)) {
    return text
  }

  return text.toLowerCase()
}

/**
 * Returns a specific number of characters from a text string starting at the position you specify
 *
 * Category: Text
 *
 * @param {*} text The text string containing the characters you want to extract.
 * @param {*} start_num The position of the first character you want to extract in text. The first character in text has start_num 1, and so on.
 * @param {*} num_chars Specifies the number of characters you want MID to return from text.
 * @returns
 */
function MID(text, start_num, num_chars) {
  if (start_num === undefined || start_num === null) {
    return value
  }

  start_num = parseNumber(start_num);
  num_chars = parseNumber(num_chars);

  if (anyIsError(start_num, num_chars) || typeof text !== 'string') {
    return num_chars
  }

  const begin = start_num - 1;
  const end = begin + num_chars;

  return text.substring(begin, end)
}

// TODO
/**
 * Converts text to number in a locale-independent manner.
 *
 * Category: Text
 *
 * @param {*} text The text to convert to a number.
 * @param {*} decimal_separator Optional. The character used to separate the integer and fractional part of the result.
 * @param {*} group_separator Optional. The character used to separate groupings of numbers, such as thousands from hundreds and millions from thousands.
 * @returns
 */
function NUMBERVALUE(text, decimal_separator, group_separator) {
  text = isDefined(text) ? text : '';

  if (typeof text === 'number') {
    return text
  }

  if (typeof text !== 'string') {
    return na
  }

  decimal_separator = typeof decimal_separator === 'undefined' ? '.' : decimal_separator;
  group_separator = typeof group_separator === 'undefined' ? ',' : group_separator;

  return Number(text.replace(decimal_separator, '.').replace(group_separator, ''))
}

// TODO
/**
 * -- Not implemented --
 */
function PRONETIC() {
  throw new Error('PRONETIC is not implemented')
}

/**
 * Capitalizes the first letter in each word of a text value.
 *
 * Category: Text
 *
 * @param {*} text Text enclosed in quotation marks, a formula that returns text, or a reference to a value containing the text you want to partially capitalize.
 * @returns
 */
function PROPER(text) {
  if (anyIsError(text)) {
    return text
  }

  if (isNaN(text) && typeof text === 'number') {
    return value
  }

  text = parseString(text);

  return text.replace(/\w\S*/g, (txt) => txt.charAt(0).toUpperCase() + txt.substr(1).toLowerCase())
}

/**
 * Replaces characters within text
 *
 * Category: Text
 *
 * @param {*} old_text Text in which you want to replace some characters.
 * @param {*} num_chars The number of characters in old_text that you want REPLACE to replace with new_text.
 * @param {*} length he number of characters in old_text that you want REPLACEB to replace with new_text.
 * @param {*} new_text he text that will replace characters in old_text.
 * @returns
 */
function REPLACE(old_text, num_chars, length, new_text) {
  num_chars = parseNumber(num_chars);
  length = parseNumber(length);

  if (anyIsError(num_chars, length) || typeof old_text !== 'string' || typeof new_text !== 'string') {
    return value
  }

  return old_text.substr(0, num_chars - 1) + new_text + old_text.substr(num_chars - 1 + length)
}

/**
 * Repeats text a given number of times.
 *
 * Category: Text
 *
 * @param {*} text The text you want to repeat.
 * @param {*} number_times A positive number specifying the number of times to repeat text.
 * @returns
 */
function REPT(text, number_times) {
  const someError = anyError(text, number_times);

  if (someError) {
    return someError
  }

  text = parseString(text);
  number_times = parseNumber(number_times);

  if (number_times instanceof Error) {
    return number_times
  }

  return new Array(number_times + 1).join(text)
}

/**
 * Returns the rightmost characters from a text value
 *
 * Category: Text
 *
 * @param {*} text The text string containing the characters you want to extract.
 * @param {*} num_chars Optional. Specifies the number of characters you want RIGHT to extract.
 * @returns
 */
function RIGHT(text, num_chars) {
  const someError = anyError(text, num_chars);

  if (someError) {
    return someError
  }

  text = parseString(text);
  num_chars = num_chars === undefined ? 1 : num_chars;
  num_chars = parseNumber(num_chars);

  if (num_chars instanceof Error) {
    return num_chars
  }

  return text.substring(text.length - num_chars)
}

/**
 * Finds one text value within another (not case-sensitive)
 *
 * Category: Text
 *
 * @param {*} find_text The text that you want to find.
 * @param {*} within_text The text in which you want to search for the value of the find_text argument.
 * @param {*} start_num Optional. The character number in the within_text argument at which you want to start searching.
 * @returns
 */
function SEARCH(find_text, within_text, start_num) {
  let foundAt;

  if (typeof find_text !== 'string' || typeof within_text !== 'string') {
    return value
  }

  start_num = start_num === undefined ? 0 : start_num;
  foundAt = within_text.toLowerCase().indexOf(find_text.toLowerCase(), start_num - 1) + 1;

  return foundAt === 0 ? value : foundAt
}

/**
 * Substitutes new text for old text in a text string.
 *
 * Category: Text
 *
 * @param {*} text The text or the reference to a value containing text for which you want to substitute characters.
 * @param {*} old_text The text you want to replace.
 * @param {*} new_text The text you want to replace old_text with.
 * @param {*} instance_num Optional. Specifies which occurrence of old_text you want to replace with new_text. If you specify instance_num, only that instance of old_text is replaced. Otherwise, every occurrence of old_text in text is changed to new_text.
 * @returns
 */
function SUBSTITUTE(text, old_text, new_text, instance_num) {
  if (arguments.length < 3) {
    return na
  }

  if (!text || !old_text) {
    return text
  } else if (instance_num === undefined) {
    return text.split(old_text).join(new_text)
  } else {
    instance_num = Math.floor(Number(instance_num));

    if (Number.isNaN(instance_num) || instance_num <= 0) {
      return value
    }

    let index = 0;
    let i = 0;

    while (index > -1 && text.indexOf(old_text, index) > -1) {
      index = text.indexOf(old_text, index + 1);
      i++;

      if (index > -1 && i === instance_num) {
        return text.substring(0, index) + new_text + text.substring(index + old_text.length)
      }
    }

    return text
  }
}

/**
 * Converts its arguments to text.
 *
 * Category: Text
 *
 * @param {*} value The value you want to test.
 * @returns
 */
function T(value) {
  if (value instanceof Error) {
    return value
  }

  return typeof value === 'string' ? value : ''
}

/**
 * Formats a number and converts it to text.
 *
 * Category: Text
 *
 * @param {*} value A numeric value that you want to be converted into text.
 * @param {*} format_text A text string that defines the formatting that you want to be applied to the supplied value.
 * @returns
 */
function TEXT(value$1, format_text) {
  if (value$1 === undefined || value$1 instanceof Error || format_text instanceof Error) {
    return na
  }

  if (format_text === undefined || format_text === null) {
    return ''
  }

  if (typeof format_text === 'number') {
    return String(format_text)
  }

  if (typeof format_text !== 'string') {
    return value
  }

  const currencySymbol = format_text.startsWith('$') ? '$' : '';
  const isPercent = format_text.endsWith('%');
  format_text = format_text.replace(/%/g, '').replace(/\$/g, '');

  // count all 0s after the decimal point
  const decimalPlaces = format_text.split('.')[1].match(/0/g).length;

  const noCommas = !format_text.includes(',');

  if (isPercent) {
    value$1 = value$1 * 100;
  }

  value$1 = FIXED(value$1, decimalPlaces, noCommas);

  if (value$1.startsWith('-')) {
    value$1 = value$1.replace('-', '');
    value$1 = '-' + currencySymbol + value$1;
  } else {
    value$1 = currencySymbol + value$1;
  }

  if (isPercent) {
    value$1 = value$1 + '%';
  }

  return value$1
}

/**
 * Combines the text from multiple ranges and/or strings.
 *
 * Category: Text
 * @param {*} delimiter A text string, either empty, or one or more characters enclosed by double quotes, or a reference to a valid text string. If a number is supplied, it will be treated as text.
 * @param {*} ignore_empty If TRUE, ignores empty values.
 * @param {*} args Text item to be joined. A text string, or array of strings, such as a range of values.
 * @returns
 */
function TEXTJOIN(delimiter, ignore_empty, ...args) {
  if (typeof ignore_empty !== 'boolean') {
    ignore_empty = parseBool(ignore_empty);
  }

  if (arguments.length < 3) {
    return na
  }

  delimiter = delimiter !== null && delimiter !== undefined ? delimiter : '';

  let flatArgs = flatten(args);
  let textToJoin = ignore_empty ? flatArgs.filter((text) => text) : flatArgs;

  if (Array.isArray(delimiter)) {
    delimiter = flatten(delimiter);

    let chunks = textToJoin.map((item) => [item]);
    let index = 0;

    for (let i = 0; i < chunks.length - 1; i++) {
      chunks[i].push(delimiter[index]);
      index++;

      if (index === delimiter.length) {
        index = 0;
      }
    }

    textToJoin = flatten(chunks);

    return textToJoin.join('')
  }

  return textToJoin.join(delimiter)
}

/**
 * Removes spaces from text.
 *
 * Category: Text
 *
 * @param {*} text The text from which you want spaces removed.
 * @returns
 */
function TRIM(text) {
  text = parseString(text);

  if (text instanceof Error) {
    return text
  }

  return text.replace(/\s+/g, ' ').trim()
}

const UNICHAR = CHAR;

const UNICODE = CODE;

/**
 * Converts text to uppercase.
 *
 * Category: Text
 *
 * @param {*} text The text you want converted to uppercase. Text can be a reference or text string.
 * @returns
 */
function UPPER(text) {
  text = parseString(text);

  if (text instanceof Error) {
    return text
  }

  return text.toUpperCase()
}

/**
 * Converts a text argument to a number.
 *
 * Category: Text
 *
 * @param {*} text The text enclosed in quotation marks or a reference to a value containing the text you want to convert.
 * @returns
 */
function VALUE(text) {
  const anyError$1 = anyError(text);

  if (anyError$1) {
    return anyError$1
  }

  if (typeof text === 'number') {
    return text
  }

  if (!isDefined(text)) {
    text = '';
  }

  if (typeof text !== 'string') {
    return value
  }

  const isPercent = /(%)$/.test(text) || /^(%)/.test(text);
  text = text.replace(/^[^0-9-]{0,3}/, '');
  text = text.replace(/[^0-9]{0,3}$/, '');
  text = text.replace(/[ ,]/g, '');

  if (text === '') {
    return 0
  }

  let output = Number(text);

  if (isNaN(output)) {
    return value
  }

  output = output || 0;

  if (isPercent) {
    output = output * 0.01;
  }

  return output
}

function isValidBinaryNumber(number) {
  return /^[01]{1,10}$/.test(number)
}

/**
 * Returns the modified Bessel function In(x).
 *
 * Category: Engineering
 *
 * @param {*} x The value at which to evaluate the function.
 * @param {*} n The order of the Bessel function. If n is not an integer, it is truncated.
 * @returns
 */
function BESSELI(x, n) {
  x = parseNumber(x);
  n = parseNumber(n);

  if (anyIsError(x, n)) {
    return value
  }

  return bessel.besseli(x, n)
}

/**
 * Returns the Bessel function Jn(x).
 *
 * Category: Engineering
 *
 * @param {*} x The value at which to evaluate the function.
 * @param {*} n The order of the Bessel function. If n is not an integer, it is truncated.
 * @returns
 */
function BESSELJ(x, n) {
  x = parseNumber(x);
  n = parseNumber(n);

  if (anyIsError(x, n)) {
    return value
  }

  return bessel.besselj(x, n)
}

/**
 * Returns the modified Bessel function Kn(x).
 *
 * Category: Engineering
 *
 * @param {*} x The value at which to evaluate the function.
 * @param {*} n The order of the function. If n is not an integer, it is truncated.
 * @returns
 */
function BESSELK(x, n) {
  x = parseNumber(x);
  n = parseNumber(n);

  if (anyIsError(x, n)) {
    return value
  }

  return bessel.besselk(x, n)
}

/**
 * Returns the Bessel function Yn(x).
 *
 * Category: Engineering
 *
 * @param {*} x The value at which to evaluate the function.
 * @param {*} n The order of the function. If n is not an integer, it is truncated.
 * @returns
 */
function BESSELY(x, n) {
  x = parseNumber(x);
  n = parseNumber(n);

  if (anyIsError(x, n)) {
    return value
  }

  return bessel.bessely(x, n)
}

/**
 * Converts a binary number to decimal.
 *
 * Category: Engineering
 *
 * @param {*} number The binary number you want to convert. Number cannot contain more than 10 characters (10 bits). The most significant bit of number is the sign bit. The remaining 9 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @returns
 */
function BIN2DEC(number) {
  // Return error if number is not binary or contains more than 10 characters (10 digits)
  if (!isValidBinaryNumber(number)) {
    return num
  }

  // Convert binary number to decimal
  const result = parseInt(number, 2);

  // Handle negative numbers
  const stringified = number.toString();

  if (stringified.length === 10 && stringified.substring(0, 1) === '1') {
    return parseInt(stringified.substring(1), 2) - 512
  } else {
    return result
  }
}

/**
 * Converts a binary number to hexadecimal.
 *
 * Category: Engineering
 *
 * @param {*} number The binary number you want to convert. Number cannot contain more than 10 characters (10 bits). The most significant bit of number is the sign bit. The remaining 9 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @param {*} places Optional. The number of characters to use. If places is omitted, BIN2HEX uses the minimum number of characters necessary. Places is useful for padding the return value with leading 0s (zeros).
 * @returns
 */
function BIN2HEX(number, places) {
  // Return error if number is not binary or contains more than 10 characters (10 digits)
  if (!isValidBinaryNumber(number)) {
    return num
  }

  // Ignore places and return a 10-character hexadecimal number if number is negative
  const stringified = number.toString();

  if (stringified.length === 10 && stringified.substring(0, 1) === '1') {
    return (1099511627264 + parseInt(stringified.substring(1), 2)).toString(16)
  }

  // Convert binary number to hexadecimal
  const result = parseInt(number, 2).toString(16);

  // Return hexadecimal number using the minimum number of characters necessary if places is undefined
  if (places === undefined) {
    return result
  } else {
    // Return error if places is nonnumeric
    if (isNaN(places)) {
      return value
    }

    // Return error if places is negative
    if (places < 0) {
      return num
    }

    // Truncate places in case it is not an integer
    places = Math.floor(places);

    // Pad return value with leading 0s (zeros) if necessary (using Underscore.string)
    return places >= result.length ? REPT('0', places - result.length) + result : num
  }
}

/**
 * Converts a binary number to octal.
 *
 * Category: Engineering
 *
 * @param {*} number The binary number you want to convert. Number cannot contain more than 10 characters (10 bits). The most significant bit of number is the sign bit. The remaining 9 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @param {*} places Optional. The number of characters to use. If places is omitted, BIN2OCT uses the minimum number of characters necessary. Places is useful for padding the return value with leading 0s (zeros).
 * @returns
 */
function BIN2OCT(number, places) {
  // Return error if number is not binary or contains more than 10 characters (10 digits)
  if (!isValidBinaryNumber(number)) {
    return num
  }

  // Ignore places and return a 10-character octal number if number is negative
  const stringified = number.toString();

  if (stringified.length === 10 && stringified.substring(0, 1) === '1') {
    return (1073741312 + parseInt(stringified.substring(1), 2)).toString(8)
  }

  // Convert binary number to octal
  const result = parseInt(number, 2).toString(8);

  // Return octal number using the minimum number of characters necessary if places is undefined
  if (places === undefined) {
    return result
  } else {
    // Return error if places is nonnumeric
    if (isNaN(places)) {
      return value
    }

    // Return error if places is negative
    if (places < 0) {
      return num
    }

    // Truncate places in case it is not an integer
    places = Math.floor(places);

    // Pad return value with leading 0s (zeros) if necessary (using Underscore.string)
    return places >= result.length ? REPT('0', places - result.length) + result : num
  }
}

/**
 * Returns a 'Bitwise And' of two numbers.
 *
 * Category: Engineering
 *
 * @param {*} number1 Must be in decimal form and greater than or equal to 0.
 * @param {*} number2 Must be in decimal form and greater than or equal to 0.
 * @returns
 */
function BITAND(number1, number2) {
  // Return error if either number is a non-numeric value
  number1 = parseNumber(number1);
  number2 = parseNumber(number2);

  if (anyIsError(number1, number2)) {
    return value
  }

  // Return error if either number is less than 0
  if (number1 < 0 || number2 < 0) {
    return num
  }

  // Return error if either number is a non-integer
  if (Math.floor(number1) !== number1 || Math.floor(number2) !== number2) {
    return num
  }

  // Return error if either number is greater than (2^48)-1
  if (number1 > 281474976710655 || number2 > 281474976710655) {
    return num
  }

  // Return bitwise AND of two numbers
  return number1 & number2
}

/**
 * Returns a value number shifted left by shift_amount bits.
 *
 * Category: Engineering
 *
 * @param {*} number Number must be an integer greater than or equal to 0.
 * @param {*} shift_amount Shift_amount must be an integer.
 * @returns
 */
function BITLSHIFT(number, shift_amount) {
  number = parseNumber(number);
  shift_amount = parseNumber(shift_amount);

  if (anyIsError(number, shift_amount)) {
    return value
  }

  // Return error if number is less than 0
  if (number < 0) {
    return num
  }

  // Return error if number is a non-integer
  if (Math.floor(number) !== number) {
    return num
  }

  // Return error if number is greater than (2^48)-1
  if (number > 281474976710655) {
    return num
  }

  // Return error if the absolute value of shift is greater than 53
  if (Math.abs(shift_amount) > 53) {
    return num
  }

  // Return number shifted by shift bits to the left or to the right if shift is negative
  return shift_amount >= 0 ? number << shift_amount : number >> -shift_amount
}

/**
 * Returns a bitwise OR of 2 numbers.
 *
 * Category: Engineering
 *
 * @param {*} number1 Must be in decimal form and greater than or equal to 0.
 * @param {*} number2 Must be in decimal form and greater than or equal to 0.
 * @returns
 */
function BITOR(number1, number2) {
  number1 = parseNumber(number1);
  number2 = parseNumber(number2);

  if (anyIsError(number1, number2)) {
    return value
  }

  // Return error if either number is less than 0
  if (number1 < 0 || number2 < 0) {
    return num
  }

  // Return error if either number is a non-integer
  if (Math.floor(number1) !== number1 || Math.floor(number2) !== number2) {
    return num
  }

  // Return error if either number is greater than (2^48)-1
  if (number1 > 281474976710655 || number2 > 281474976710655) {
    return num
  }

  // Return bitwise OR of two numbers
  return number1 | number2
}

/**
 * Returns a value number shifted right by shift_amount bits.
 *
 * Category: Engineering
 *
 * @param {*} number Must be an integer greater than or equal to 0.
 * @param {*} shift_amount Must be an integer.
 * @returns
 */
function BITRSHIFT(number, shift_amount) {
  number = parseNumber(number);
  shift_amount = parseNumber(shift_amount);

  if (anyIsError(number, shift_amount)) {
    return value
  }

  // Return error if number is less than 0
  if (number < 0) {
    return num
  }

  // Return error if number is a non-integer
  if (Math.floor(number) !== number) {
    return num
  }

  // Return error if number is greater than (2^48)-1
  if (number > 281474976710655) {
    return num
  }

  // Return error if the absolute value of shift is greater than 53
  if (Math.abs(shift_amount) > 53) {
    return num
  }

  // Return number shifted by shift bits to the right or to the left if shift is negative
  return shift_amount >= 0 ? number >> shift_amount : number << -shift_amount
}

/**
 * Returns a bitwise 'Exclusive Or' of two numbers.
 *
 * Category: Engineering
 *
 * @param {*} number1 Must be greater than or equal to 0.
 * @param {*} number2 Must be greater than or equal to 0.
 * @returns
 */
function BITXOR(number1, number2) {
  number1 = parseNumber(number1);
  number2 = parseNumber(number2);

  if (anyIsError(number1, number2)) {
    return value
  }

  // Return error if either number is less than 0
  if (number1 < 0 || number2 < 0) {
    return num
  }

  // Return error if either number is a non-integer
  if (Math.floor(number1) !== number1 || Math.floor(number2) !== number2) {
    return num
  }

  // Return error if either number is greater than (2^48)-1
  if (number1 > 281474976710655 || number2 > 281474976710655) {
    return num
  }

  // Return bitwise XOR of two numbers
  return number1 ^ number2
}

/**
 * Converts real and imaginary coefficients into a complex number.
 *
 * Category: Engineering
 *
 * @param {*} real_num The real coefficient of the complex number.
 * @param {*} i_num The imaginary coefficient of the complex number.
 * @param {*} suffix Optional. The suffix for the imaginary component of the complex number. If omitted, suffix is assumed to be "i".
 * @returns
 */
function COMPLEX(real_num, i_num, suffix) {
  real_num = parseNumber(real_num);
  i_num = parseNumber(i_num);

  if (anyIsError(real_num, i_num)) {
    return real_num
  }

  // Set suffix
  suffix = suffix === undefined ? 'i' : suffix;

  // Return error if suffix is neither "i" nor "j"
  if (suffix !== 'i' && suffix !== 'j') {
    return value
  }

  // Return complex number
  if (real_num === 0 && i_num === 0) {
    return 0
  } else if (real_num === 0) {
    return i_num === 1 ? suffix : i_num.toString() + suffix
  } else if (i_num === 0) {
    return real_num.toString()
  } else {
    const sign = i_num > 0 ? '+' : '';
    return real_num.toString() + sign + (i_num === 1 ? suffix : i_num.toString() + suffix)
  }
}

/**
 * Converts a number from one measurement system to another.
 *
 * Category: Engineering
 *
 * @param {*} number is the value in from_units to convert.
 * @param {*} from_unit is the units for number.
 * @param {*} to_unit is the units for the result. CONVERT accepts the following text values (in quotation marks) for from_unit and to_unit.
 * @returns
 */
function CONVERT(number, from_unit, to_unit) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  // List of units supported by CONVERT and units defined by the International System of Units
  // [Name, Symbol, Alternate symbols, Quantity, ISU, CONVERT, Conversion ratio]
  const units = [
    ['a.u. of action', '?', null, 'action', false, false, 1.05457168181818e-34],
    ['a.u. of charge', 'e', null, 'electric_charge', false, false, 1.60217653141414e-19],
    ['a.u. of energy', 'Eh', null, 'energy', false, false, 4.35974417757576e-18],
    ['a.u. of length', 'a?', null, 'length', false, false, 5.29177210818182e-11],
    ['a.u. of mass', 'm?', null, 'mass', false, false, 9.10938261616162e-31],
    ['a.u. of time', '?/Eh', null, 'time', false, false, 2.41888432650516e-17],
    ['admiralty knot', 'admkn', null, 'speed', false, true, 0.514773333],
    ['ampere', 'A', null, 'electric_current', true, false, 1],
    ['ampere per meter', 'A/m', null, 'magnetic_field_intensity', true, false, 1],
    ['ångström', 'Å', ['ang'], 'length', false, true, 1e-10],
    ['are', 'ar', null, 'area', false, true, 100],
    ['astronomical unit', 'ua', null, 'length', false, false, 1.49597870691667e-11],
    ['bar', 'bar', null, 'pressure', false, false, 100000],
    ['barn', 'b', null, 'area', false, false, 1e-28],
    ['becquerel', 'Bq', null, 'radioactivity', true, false, 1],
    ['bit', 'bit', ['b'], 'information', false, true, 1],
    ['btu', 'BTU', ['btu'], 'energy', false, true, 1055.05585262],
    ['byte', 'byte', null, 'information', false, true, 8],
    ['candela', 'cd', null, 'luminous_intensity', true, false, 1],
    ['candela per square metre', 'cd/m?', null, 'luminance', true, false, 1],
    ['coulomb', 'C', null, 'electric_charge', true, false, 1],
    ['cubic ångström', 'ang3', ['ang^3'], 'volume', false, true, 1e-30],
    ['cubic foot', 'ft3', ['ft^3'], 'volume', false, true, 0.028316846592],
    ['cubic inch', 'in3', ['in^3'], 'volume', false, true, 0.000016387064],
    ['cubic light-year', 'ly3', ['ly^3'], 'volume', false, true, 8.46786664623715e-47],
    ['cubic metre', 'm?', null, 'volume', true, true, 1],
    ['cubic mile', 'mi3', ['mi^3'], 'volume', false, true, 4168181825.44058],
    ['cubic nautical mile', 'Nmi3', ['Nmi^3'], 'volume', false, true, 6352182208],
    ['cubic Pica', 'Pica3', ['Picapt3', 'Pica^3', 'Picapt^3'], 'volume', false, true, 7.58660370370369e-8],
    ['cubic yard', 'yd3', ['yd^3'], 'volume', false, true, 0.764554857984],
    ['cup', 'cup', null, 'volume', false, true, 0.0002365882365],
    ['dalton', 'Da', ['u'], 'mass', false, false, 1.66053886282828e-27],
    ['day', 'd', ['day'], 'time', false, true, 86400],
    ['degree', '°', null, 'angle', false, false, 0.0174532925199433],
    ['degrees Rankine', 'Rank', null, 'temperature', false, true, 0.555555555555556],
    ['dyne', 'dyn', ['dy'], 'force', false, true, 0.00001],
    ['electronvolt', 'eV', ['ev'], 'energy', false, true, 1.60217656514141],
    ['ell', 'ell', null, 'length', false, true, 1.143],
    ['erg', 'erg', ['e'], 'energy', false, true, 1e-7],
    ['farad', 'F', null, 'electric_capacitance', true, false, 1],
    ['fluid ounce', 'oz', null, 'volume', false, true, 0.0000295735295625],
    ['foot', 'ft', null, 'length', false, true, 0.3048],
    ['foot-pound', 'flb', null, 'energy', false, true, 1.3558179483314],
    ['gal', 'Gal', null, 'acceleration', false, false, 0.01],
    ['gallon', 'gal', null, 'volume', false, true, 0.003785411784],
    ['gauss', 'G', ['ga'], 'magnetic_flux_density', false, true, 1],
    ['grain', 'grain', null, 'mass', false, true, 0.0000647989],
    ['gram', 'g', null, 'mass', false, true, 0.001],
    ['gray', 'Gy', null, 'absorbed_dose', true, false, 1],
    ['gross registered ton', 'GRT', ['regton'], 'volume', false, true, 2.8316846592],
    ['hectare', 'ha', null, 'area', false, true, 10000],
    ['henry', 'H', null, 'inductance', true, false, 1],
    ['hertz', 'Hz', null, 'frequency', true, false, 1],
    ['horsepower', 'HP', ['h'], 'power', false, true, 745.69987158227],
    ['horsepower-hour', 'HPh', ['hh', 'hph'], 'energy', false, true, 2684519.538],
    ['hour', 'h', ['hr'], 'time', false, true, 3600],
    ['imperial gallon (U.K.)', 'uk_gal', null, 'volume', false, true, 0.00454609],
    ['imperial hundredweight', 'lcwt', ['uk_cwt', 'hweight'], 'mass', false, true, 50.802345],
    ['imperial quart (U.K)', 'uk_qt', null, 'volume', false, true, 0.0011365225],
    ['imperial ton', 'brton', ['uk_ton', 'LTON'], 'mass', false, true, 1016.046909],
    ['inch', 'in', null, 'length', false, true, 0.0254],
    ['international acre', 'uk_acre', null, 'area', false, true, 4046.8564224],
    ['IT calorie', 'cal', null, 'energy', false, true, 4.1868],
    ['joule', 'J', null, 'energy', true, true, 1],
    ['katal', 'kat', null, 'catalytic_activity', true, false, 1],
    ['kelvin', 'K', ['kel'], 'temperature', true, true, 1],
    ['kilogram', 'kg', null, 'mass', true, true, 1],
    ['knot', 'kn', null, 'speed', false, true, 0.514444444444444],
    ['light-year', 'ly', null, 'length', false, true, 9460730472580800],
    ['litre', 'L', ['l', 'lt'], 'volume', false, true, 0.001],
    ['lumen', 'lm', null, 'luminous_flux', true, false, 1],
    ['lux', 'lx', null, 'illuminance', true, false, 1],
    ['maxwell', 'Mx', null, 'magnetic_flux', false, false, 1e-18],
    ['measurement ton', 'MTON', null, 'volume', false, true, 1.13267386368],
    ['meter per hour', 'm/h', ['m/hr'], 'speed', false, true, 0.00027777777777778],
    ['meter per second', 'm/s', ['m/sec'], 'speed', true, true, 1],
    ['meter per second squared', 'm?s??', null, 'acceleration', true, false, 1],
    ['parsec', 'pc', ['parsec'], 'length', false, true, 30856775814671900],
    ['meter squared per second', 'm?/s', null, 'kinematic_viscosity', true, false, 1],
    ['metre', 'm', null, 'length', true, true, 1],
    ['miles per hour', 'mph', null, 'speed', false, true, 0.44704],
    ['millimetre of mercury', 'mmHg', null, 'pressure', false, false, 133.322],
    ['minute', '?', null, 'angle', false, false, 0.000290888208665722],
    ['minute', 'min', ['mn'], 'time', false, true, 60],
    ['modern teaspoon', 'tspm', null, 'volume', false, true, 0.000005],
    ['mole', 'mol', null, 'amount_of_substance', true, false, 1],
    ['morgen', 'Morgen', null, 'area', false, true, 2500],
    ['n.u. of action', '?', null, 'action', false, false, 1.05457168181818e-34],
    ['n.u. of mass', 'm?', null, 'mass', false, false, 9.10938261616162e-31],
    ['n.u. of speed', 'c?', null, 'speed', false, false, 299792458],
    ['n.u. of time', '?/(me?c??)', null, 'time', false, false, 1.28808866778687e-21],
    ['nautical mile', 'M', ['Nmi'], 'length', false, true, 1852],
    ['newton', 'N', null, 'force', true, true, 1],
    ['œrsted', 'Oe ', null, 'magnetic_field_intensity', false, false, 79.5774715459477],
    ['ohm', 'Ω', null, 'electric_resistance', true, false, 1],
    ['ounce mass', 'ozm', null, 'mass', false, true, 0.028349523125],
    ['pascal', 'Pa', null, 'pressure', true, false, 1],
    ['pascal second', 'Pa?s', null, 'dynamic_viscosity', true, false, 1],
    ['pferdestärke', 'PS', null, 'power', false, true, 735.49875],
    ['phot', 'ph', null, 'illuminance', false, false, 0.0001],
    ['pica (1/6 inch)', 'pica', null, 'length', false, true, 0.00035277777777778],
    ['pica (1/72 inch)', 'Pica', ['Picapt'], 'length', false, true, 0.00423333333333333],
    ['poise', 'P', null, 'dynamic_viscosity', false, false, 0.1],
    ['pond', 'pond', null, 'force', false, true, 0.00980665],
    ['pound force', 'lbf', null, 'force', false, true, 4.4482216152605],
    ['pound mass', 'lbm', null, 'mass', false, true, 0.45359237],
    ['quart', 'qt', null, 'volume', false, true, 0.000946352946],
    ['radian', 'rad', null, 'angle', true, false, 1],
    ['second', '?', null, 'angle', false, false, 0.00000484813681109536],
    ['second', 's', ['sec'], 'time', true, true, 1],
    ['short hundredweight', 'cwt', ['shweight'], 'mass', false, true, 45.359237],
    ['siemens', 'S', null, 'electrical_conductance', true, false, 1],
    ['sievert', 'Sv', null, 'equivalent_dose', true, false, 1],
    ['slug', 'sg', null, 'mass', false, true, 14.59390294],
    ['square ångström', 'ang2', ['ang^2'], 'area', false, true, 1e-20],
    ['square foot', 'ft2', ['ft^2'], 'area', false, true, 0.09290304],
    ['square inch', 'in2', ['in^2'], 'area', false, true, 0.00064516],
    ['square light-year', 'ly2', ['ly^2'], 'area', false, true, 8.95054210748189e31],
    ['square meter', 'm?', null, 'area', true, true, 1],
    ['square mile', 'mi2', ['mi^2'], 'area', false, true, 2589988.110336],
    ['square nautical mile', 'Nmi2', ['Nmi^2'], 'area', false, true, 3429904],
    ['square Pica', 'Pica2', ['Picapt2', 'Pica^2', 'Picapt^2'], 'area', false, true, 0.00001792111111111],
    ['square yard', 'yd2', ['yd^2'], 'area', false, true, 0.83612736],
    ['statute mile', 'mi', null, 'length', false, true, 1609.344],
    ['steradian', 'sr', null, 'solid_angle', true, false, 1],
    ['stilb', 'sb', null, 'luminance', false, false, 0.0001],
    ['stokes', 'St', null, 'kinematic_viscosity', false, false, 0.0001],
    ['stone', 'stone', null, 'mass', false, true, 6.35029318],
    ['tablespoon', 'tbs', null, 'volume', false, true, 0.0000147868],
    ['teaspoon', 'tsp', null, 'volume', false, true, 0.00000492892],
    ['tesla', 'T', null, 'magnetic_flux_density', true, true, 1],
    ['thermodynamic calorie', 'c', null, 'energy', false, true, 4.184],
    ['ton', 'ton', null, 'mass', false, true, 907.18474],
    ['tonne', 't', null, 'mass', false, false, 1000],
    ['U.K. pint', 'uk_pt', null, 'volume', false, true, 0.00056826125],
    ['U.S. bushel', 'bushel', null, 'volume', false, true, 0.03523907],
    ['U.S. oil barrel', 'barrel', null, 'volume', false, true, 0.158987295],
    ['U.S. pint', 'pt', ['us_pt'], 'volume', false, true, 0.000473176473],
    ['U.S. survey mile', 'survey_mi', null, 'length', false, true, 1609.347219],
    ['U.S. survey/statute acre', 'us_acre', null, 'area', false, true, 4046.87261],
    ['volt', 'V', null, 'voltage', true, false, 1],
    ['watt', 'W', null, 'power', true, true, 1],
    ['watt-hour', 'Wh', ['wh'], 'energy', false, true, 3600],
    ['weber', 'Wb', null, 'magnetic_flux', true, false, 1],
    ['yard', 'yd', null, 'length', false, true, 0.9144],
    ['year', 'yr', null, 'time', false, true, 31557600]
  ];

  // Binary prefixes
  // [Name, Prefix power of 2 value, Previx value, Abbreviation, Derived from]
  const binary_prefixes = {
    Yi: ['yobi', 80, 1208925819614629174706176, 'Yi', 'yotta'],
    Zi: ['zebi', 70, 1180591620717411303424, 'Zi', 'zetta'],
    Ei: ['exbi', 60, 1152921504606846976, 'Ei', 'exa'],
    Pi: ['pebi', 50, 1125899906842624, 'Pi', 'peta'],
    Ti: ['tebi', 40, 1099511627776, 'Ti', 'tera'],
    Gi: ['gibi', 30, 1073741824, 'Gi', 'giga'],
    Mi: ['mebi', 20, 1048576, 'Mi', 'mega'],
    ki: ['kibi', 10, 1024, 'ki', 'kilo']
  };

  // Unit prefixes
  // [Name, Multiplier, Abbreviation]
  const unit_prefixes = {
    Y: ['yotta', 1e24, 'Y'],
    Z: ['zetta', 1e21, 'Z'],
    E: ['exa', 1e18, 'E'],
    P: ['peta', 1e15, 'P'],
    T: ['tera', 1e12, 'T'],
    G: ['giga', 1e9, 'G'],
    M: ['mega', 1e6, 'M'],
    k: ['kilo', 1e3, 'k'],
    h: ['hecto', 1e2, 'h'],
    e: ['dekao', 1e1, 'e'],
    d: ['deci', 1e-1, 'd'],
    c: ['centi', 1e-2, 'c'],
    m: ['milli', 1e-3, 'm'],
    u: ['micro', 1e-6, 'u'],
    n: ['nano', 1e-9, 'n'],
    p: ['pico', 1e-12, 'p'],
    f: ['femto', 1e-15, 'f'],
    a: ['atto', 1e-18, 'a'],
    z: ['zepto', 1e-21, 'z'],
    y: ['yocto', 1e-24, 'y']
  };

  // Initialize units and multipliers
  let from = null;
  let to = null;
  let base_from_unit = from_unit;
  let base_to_unit = to_unit;
  let from_multiplier = 1;
  let to_multiplier = 1;
  let alt;

  // Lookup from and to units
  for (let i = 0; i < units.length; i++) {
    alt = units[i][2] === null ? [] : units[i][2];

    if (units[i][1] === base_from_unit || alt.indexOf(base_from_unit) >= 0) {
      from = units[i];
    }

    if (units[i][1] === base_to_unit || alt.indexOf(base_to_unit) >= 0) {
      to = units[i];
    }
  }

  // Lookup from prefix
  if (from === null) {
    const from_binary_prefix = binary_prefixes[from_unit.substring(0, 2)];
    let from_unit_prefix = unit_prefixes[from_unit.substring(0, 1)];

    // Handle dekao unit prefix (only unit prefix with two characters)
    if (from_unit.substring(0, 2) === 'da') {
      from_unit_prefix = ['dekao', 1e1, 'da'];
    }

    // Handle binary prefixes first (so that 'Yi' is processed before 'Y')
    if (from_binary_prefix) {
      from_multiplier = from_binary_prefix[2];
      base_from_unit = from_unit.substring(2);
    } else if (from_unit_prefix) {
      from_multiplier = from_unit_prefix[1];
      base_from_unit = from_unit.substring(from_unit_prefix[2].length);
    }

    // Lookup from unit
    for (let j = 0; j < units.length; j++) {
      alt = units[j][2] === null ? [] : units[j][2];

      if (units[j][1] === base_from_unit || alt.indexOf(base_from_unit) >= 0) {
        from = units[j];
      }
    }
  }

  // Lookup to prefix
  if (to === null) {
    const to_binary_prefix = binary_prefixes[to_unit.substring(0, 2)];
    let to_unit_prefix = unit_prefixes[to_unit.substring(0, 1)];

    // Handle dekao unit prefix (only unit prefix with two characters)
    if (to_unit.substring(0, 2) === 'da') {
      to_unit_prefix = ['dekao', 1e1, 'da'];
    }

    // Handle binary prefixes first (so that 'Yi' is processed before 'Y')
    if (to_binary_prefix) {
      to_multiplier = to_binary_prefix[2];
      base_to_unit = to_unit.substring(2);
    } else if (to_unit_prefix) {
      to_multiplier = to_unit_prefix[1];
      base_to_unit = to_unit.substring(to_unit_prefix[2].length);
    }

    // Lookup to unit
    for (let k = 0; k < units.length; k++) {
      alt = units[k][2] === null ? [] : units[k][2];

      if (units[k][1] === base_to_unit || alt.indexOf(base_to_unit) >= 0) {
        to = units[k];
      }
    }
  }

  // Return error if a unit does not exist
  if (from === null || to === null) {
    return na
  }

  // Return error if units represent different quantities
  if (from[3] !== to[3]) {
    return na
  }

  // Return converted number
  return (number * from[6] * from_multiplier) / (to[6] * to_multiplier)
}

/**
 * Converts a decimal number to binary.
 *
 * Category: Engineering
 *
 * @param {*} number The decimal integer you want to convert. If number is negative, valid place values are ignored and DEC2BIN returns a 10-character (10-bit) binary number in which the most significant bit is the sign bit. The remaining 9 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @param {*} places Optional. The number of characters to use. If places is omitted, DEC2BIN uses the minimum number of characters necessary. Places is useful for padding the return value with leading 0s (zeros).
 * @returns
 */
function DEC2BIN(number, places) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  // Return error if number is not decimal, is lower than -512, or is greater than 511
  if (!/^-?[0-9]{1,3}$/.test(number) || number < -512 || number > 511) {
    return num
  }

  // Ignore places and return a 10-character binary number if number is negative
  if (number < 0) {
    return '1' + REPT('0', 9 - (512 + number).toString(2).length) + (512 + number).toString(2)
  }

  // Convert decimal number to binary
  const result = parseInt(number, 10).toString(2);

  // Return binary number using the minimum number of characters necessary if places is undefined
  if (typeof places === 'undefined') {
    return result
  } else {
    // Return error if places is nonnumeric
    if (isNaN(places)) {
      return value
    }

    // Return error if places is negative
    if (places < 0) {
      return num
    }

    // Truncate places in case it is not an integer
    places = Math.floor(places);

    // Pad return value with leading 0s (zeros) if necessary (using Underscore.string)
    return places >= result.length ? REPT('0', places - result.length) + result : num
  }
}

/**
 * Converts a decimal number to hexadecimal.
 *
 * Category: Engineering
 *
 * @param {*} number The decimal integer you want to convert. If number is negative, places is ignored and DEC2HEX returns a 10-character (40-bit) hexadecimal number in which the most significant bit is the sign bit. The remaining 39 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @param {*} places Optional. The number of characters to use. If places is omitted, DEC2HEX uses the minimum number of characters necessary. Places is useful for padding the return value with leading 0s (zeros).
 * @returns
 */
function DEC2HEX(number, places) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  // Return error if number is not decimal, is lower than -549755813888, or is greater than 549755813887
  if (!/^-?[0-9]{1,12}$/.test(number) || number < -549755813888 || number > 549755813887) {
    return num
  }

  // Ignore places and return a 10-character hexadecimal number if number is negative
  if (number < 0) {
    return (1099511627776 + number).toString(16)
  }

  // Convert decimal number to hexadecimal
  const result = parseInt(number, 10).toString(16);

  // Return hexadecimal number using the minimum number of characters necessary if places is undefined
  if (typeof places === 'undefined') {
    return result
  } else {
    // Return error if places is nonnumeric
    if (isNaN(places)) {
      return value
    }

    // Return error if places is negative
    if (places < 0) {
      return num
    }

    // Truncate places in case it is not an integer
    places = Math.floor(places);

    // Pad return value with leading 0s (zeros) if necessary (using Underscore.string)
    return places >= result.length ? REPT('0', places - result.length) + result : num
  }
}

/**
 * Converts a decimal number to octal.
 *
 * Category: Engineering
 *
 * @param {*} number The decimal integer you want to convert. If number is negative, places is ignored and DEC2OCT returns a 10-character (30-bit) octal number in which the most significant bit is the sign bit. The remaining 29 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @param {*} places Optional. The number of characters to use. If places is omitted, DEC2OCT uses the minimum number of characters necessary. Places is useful for padding the return value with leading 0s (zeros).
 * @returns
 */
function DEC2OCT(number, places) {
  number = parseNumber(number);

  if (number instanceof Error) {
    return number
  }

  // Return error if number is not decimal, is lower than -549755813888, or is greater than 549755813887
  if (!/^-?[0-9]{1,9}$/.test(number) || number < -536870912 || number > 536870911) {
    return num
  }

  // Ignore places and return a 10-character octal number if number is negative
  if (number < 0) {
    return (1073741824 + number).toString(8)
  }

  // Convert decimal number to octal
  const result = parseInt(number, 10).toString(8);

  // Return octal number using the minimum number of characters necessary if places is undefined
  if (typeof places === 'undefined') {
    return result
  } else {
    // Return error if places is nonnumeric
    if (isNaN(places)) {
      return value
    }

    // Return error if places is negative
    if (places < 0) {
      return num
    }

    // Truncate places in case it is not an integer
    places = Math.floor(places);

    // Pad return value with leading 0s (zeros) if necessary (using Underscore.string)
    return places >= result.length ? REPT('0', places - result.length) + result : num
  }
}

/**
 * Tests whether two values are equal.
 *
 * Category: Engineering
 *
 * @param {*} number1 The first number.
 * @param {*} number2 Optional. The second number. If omitted, number2 is assumed to be zero.
 * @returns
 */
function DELTA(number1, number2) {
  // Set number2 to zero if undefined
  number2 = number2 === undefined ? 0 : number2;
  number1 = parseNumber(number1);
  number2 = parseNumber(number2);

  if (anyIsError(number1, number2)) {
    return value
  }

  // Return delta
  return number1 === number2 ? 1 : 0
}

// TODO: why is upper_bound not used ? The excel documentation has no examples with upper_bound
/**
 * Returns the error function.
 *
 * Category: Engineering
 *
 * @param {*} lower_limit The lower bound for integrating ERF.
 * @param {*} upper_limit Optional. The upper bound for integrating ERF. If omitted, ERF integrates between zero and lower_limit.
 * @returns
 */
function ERF(lower_limit, upper_limit) {
  // Set number2 to zero if undefined
  upper_limit = upper_limit === undefined ? 0 : upper_limit;

  lower_limit = parseNumber(lower_limit);
  upper_limit = parseNumber(upper_limit);

  if (anyIsError(lower_limit, upper_limit)) {
    return value
  }

  return jstat.erf(lower_limit)
}

// TODO

/**
 * -- Not implemented --
 *
 * Returns the error function.
 *
 * Category: Engineering
 *
 * @param {*} x The lower bound for integrating ERF.PRECISE.
 * @returns
 */
ERF.PRECISE = () => {
  throw new Error('ERF.PRECISE is not implemented')
};

/**
 * Returns the complementary error function.
 *
 * Category: Engineering
 *
 * @param {*} x The lower bound for integrating ERFC.
 * @returns
 */
function ERFC(x) {
  // Return error if x is not a number
  if (isNaN(x)) {
    return value
  }

  return jstat.erfc(x)
}

// TODO

/**
 * -- Not implemented --
 *
 * Returns the complementary ERF function integrated between x and infinity.
 *
 * Category: Engineering
 *
 * @param {*} x The lower bound for integrating ERFC.PRECISE.
 * @returns
 */
ERFC.PRECISE = () => {
  throw new Error('ERFC.PRECISE is not implemented')
};

/**
 * Tests whether a number is greater than a threshold value.
 *
 * Category: Engineering
 *
 * @param {*} number The value to test against step.
 * @param {*} step Optional. The threshold value. If you omit a value for step, GESTEP uses zero.
 * @returns
 */
function GESTEP(number, step) {
  step = step || 0;
  number = parseNumber(number);

  if (anyIsError(step, number)) {
    return number
  }

  // Return delta
  return number >= step ? 1 : 0
}

/**
 * Converts a hexadecimal number to binary.
 *
 * Category: Engineering
 *
 * @param {*} number The hexadecimal number you want to convert. Number cannot contain more than 10 characters. The most significant bit of number is the sign bit (40th bit from the right). The remaining 9 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @param {*} places Optional. The number of characters to use. If places is omitted, HEX2BIN uses the minimum number of characters necessary. Places is useful for padding the return value with leading 0s (zeros).
 * @returns
 */
function HEX2BIN(number, places) {
  // Return error if number is not hexadecimal or contains more than ten characters (10 digits)
  if (!/^[0-9A-Fa-f]{1,10}$/.test(number)) {
    return num
  }

  // Check if number is negative
  const negative = !!(number.length === 10 && number.substring(0, 1).toLowerCase() === 'f');

  // Convert hexadecimal number to decimal
  const decimal = negative ? parseInt(number, 16) - 1099511627776 : parseInt(number, 16);

  // Return error if number is lower than -512 or greater than 511
  if (decimal < -512 || decimal > 511) {
    return num
  }

  // Ignore places and return a 10-character binary number if number is negative
  if (negative) {
    return '1' + REPT('0', 9 - (512 + decimal).toString(2).length) + (512 + decimal).toString(2)
  }

  // Convert decimal number to binary
  const result = decimal.toString(2);

  // Return binary number using the minimum number of characters necessary if places is undefined
  if (places === undefined) {
    return result
  } else {
    // Return error if places is nonnumeric
    if (isNaN(places)) {
      return value
    }

    // Return error if places is negative
    if (places < 0) {
      return num
    }

    // Truncate places in case it is not an integer
    places = Math.floor(places);

    // Pad return value with leading 0s (zeros) if necessary (using Underscore.string)
    return places >= result.length ? REPT('0', places - result.length) + result : num
  }
}

/**
 * Converts a hexadecimal number to decimal.
 *
 * Category: Engineering
 *
 * @param {*} number The hexadecimal number you want to convert. Number cannot contain more than 10 characters (40 bits). The most significant bit of number is the sign bit. The remaining 39 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @returns
 */
function HEX2DEC(number) {
  // Return error if number is not hexadecimal or contains more than ten characters (10 digits)
  if (!/^[0-9A-Fa-f]{1,10}$/.test(number)) {
    return num
  }

  // Convert hexadecimal number to decimal
  const decimal = parseInt(number, 16);

  // Return decimal number
  return decimal >= 549755813888 ? decimal - 1099511627776 : decimal
}

/**
 * Converts a hexadecimal number to octal.
 *
 * Category: Engineering
 *
 * @param {*} number The hexadecimal number you want to convert. Number cannot contain more than 10 characters. The most significant bit of number is the sign bit. The remaining 39 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @param {*} places Optional. The number of characters to use. If places is omitted, HEX2OCT uses the minimum number of characters necessary. Places is useful for padding the return value with leading 0s (zeros).
 * @returns
 */
function HEX2OCT(number, places) {
  // Return error if number is not hexadecimal or contains more than ten characters (10 digits)
  if (!/^[0-9A-Fa-f]{1,10}$/.test(number)) {
    return num
  }

  // Convert hexadecimal number to decimal
  const decimal = parseInt(number, 16);

  // Return error if number is positive and greater than 0x1fffffff (536870911)
  if (decimal > 536870911 && decimal < 1098974756864) {
    return num
  }

  // Ignore places and return a 10-character octal number if number is negative
  if (decimal >= 1098974756864) {
    return (decimal - 1098437885952).toString(8)
  }

  // Convert decimal number to octal
  const result = decimal.toString(8);

  // Return octal number using the minimum number of characters necessary if places is undefined
  if (places === undefined) {
    return result
  } else {
    // Return error if places is nonnumeric
    if (isNaN(places)) {
      return value
    }

    // Return error if places is negative
    if (places < 0) {
      return num
    }

    // Truncate places in case it is not an integer
    places = Math.floor(places);

    // Pad return value with leading 0s (zeros) if necessary (using Underscore.string)
    return places >= result.length ? REPT('0', places - result.length) + result : num
  }
}

/**
 * Returns the absolute value (modulus) of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the absolute value.
 * @returns
 */
function IMABS(inumber) {
  // Lookup real and imaginary coefficients using exports.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  // Return error if either coefficient is not a number
  if (anyIsError(x, y)) {
    return value
  }

  // Return absolute value of complex number
  return Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2))
}

/**
 * Returns the imaginary coefficient of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the imaginary coefficient.
 * @returns
 */
function IMAGINARY(inumber) {
  if (inumber === undefined || inumber === true || inumber === false) {
    return value
  }

  // Return 0 if inumber is equal to 0
  if (inumber === 0 || inumber === '0') {
    return 0
  }

  // Handle special cases
  if (['i', 'j'].indexOf(inumber) >= 0) {
    return 1
  }

  // Force string type
  inumber = inumber + '';

  // Normalize imaginary coefficient
  inumber = inumber.replace('+i', '+1i').replace('-i', '-1i').replace('+j', '+1j').replace('-j', '-1j');

  // Lookup sign
  let plus = inumber.indexOf('+');
  let minus = inumber.indexOf('-');

  if (plus === 0) {
    plus = inumber.indexOf('+', 1);
  }

  if (minus === 0) {
    minus = inumber.indexOf('-', 1);
  }

  // Lookup imaginary unit
  const last = inumber.substring(inumber.length - 1, inumber.length);
  const unit = last === 'i' || last === 'j';

  if (plus >= 0 || minus >= 0) {
    // Return error if imaginary unit is neither i nor j
    if (!unit) {
      return num
    }

    // Return imaginary coefficient of complex number
    if (plus >= 0) {
      return isNaN(inumber.substring(0, plus)) || isNaN(inumber.substring(plus + 1, inumber.length - 1))
        ? num
        : Number(inumber.substring(plus + 1, inumber.length - 1))
    } else {
      return isNaN(inumber.substring(0, minus)) || isNaN(inumber.substring(minus + 1, inumber.length - 1))
        ? num
        : -Number(inumber.substring(minus + 1, inumber.length - 1))
    }
  } else {
    if (unit) {
      return isNaN(inumber.substring(0, inumber.length - 1)) ? num : inumber.substring(0, inumber.length - 1)
    } else {
      return isNaN(inumber) ? num : 0
    }
  }
}

/**
 * Returns the argument theta, an angle expressed in radians.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the argument .
 * @returns
 */
function IMARGUMENT(inumber) {
  // Lookup real and imaginary coefficients using exports.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  // Return error if either coefficient is not a number
  if (anyIsError(x, y)) {
    return value
  }

  // Return error if inumber is equal to zero
  if (x === 0 && y === 0) {
    return div0
  }

  // Return PI/2 if x is equal to zero and y is positive
  if (x === 0 && y > 0) {
    return Math.PI / 2
  }

  // Return -PI/2 if x is equal to zero and y is negative
  if (x === 0 && y < 0) {
    return -Math.PI / 2
  }

  // Return zero if x is negative and y is equal to zero
  if (y === 0 && x > 0) {
    return 0
  }

  // Return zero if x is negative and y is equal to zero
  if (y === 0 && x < 0) {
    return -Math.PI
  }

  // Return argument of complex number
  if (x > 0) {
    return Math.atan(y / x)
  } else if (x < 0 && y >= 0) {
    return Math.atan(y / x) + Math.PI
  } else {
    return Math.atan(y / x) - Math.PI
  }
}

/**
 * Returns the complex conjugate of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the conjugate.
 * @returns
 */
function IMCONJUGATE(inumber) {
  // Lookup real and imaginary coefficients using exports.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Return conjugate of complex number
  return y !== 0 ? COMPLEX(x, -y, unit) : inumber
}

/**
 * Returns the cosine of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the cosine.
 * @returns
 */
function IMCOS(inumber) {
  // Lookup real and imaginary coefficients using exports.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Return cosine of complex number
  return COMPLEX(
    (Math.cos(x) * (Math.exp(y) + Math.exp(-y))) / 2,
    (-Math.sin(x) * (Math.exp(y) - Math.exp(-y))) / 2,
    unit
  )
}

/**
 * Returns the hyperbolic cosine of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the hyperbolic cosine.
 * @returns
 */
function IMCOSH(inumber) {
  // Lookup real and imaginary coefficients using exports.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Return hyperbolic cosine of complex number
  return COMPLEX(
    (Math.cos(y) * (Math.exp(x) + Math.exp(-x))) / 2,
    (Math.sin(y) * (Math.exp(x) - Math.exp(-x))) / 2,
    unit
  )
}

/**
 * Returns the cotangent of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the cotangent.
 * @returns
 */
function IMCOT(inumber) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Return cotangent of complex number
  return IMDIV(IMCOS(inumber), IMSIN(inumber))
}

/**
 * Returns the quotient of two complex numbers.
 *
 * Category: Engineering
 *
 * @param {*} inumber1 The complex numerator or dividend.
 * @param {*} inumber2 The complex denominator or divisor.
 * @returns
 */
function IMDIV(inumber1, inumber2) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const a = IMREAL(inumber1);
  const b = IMAGINARY(inumber1);
  const c = IMREAL(inumber2);
  const d = IMAGINARY(inumber2);

  if (anyIsError(a, b, c, d)) {
    return value
  }

  // Lookup imaginary unit
  const unit1 = inumber1.substring(inumber1.length - 1);
  const unit2 = inumber2.substring(inumber2.length - 1);
  let unit = 'i';

  if (unit1 === 'j') {
    unit = 'j';
  } else if (unit2 === 'j') {
    unit = 'j';
  }

  // Return error if inumber2 is null
  if (c === 0 && d === 0) {
    return num
  }

  // Return exponential of complex number
  const den = c * c + d * d;
  return COMPLEX((a * c + b * d) / den, (b * c - a * d) / den, unit)
}

/**
 * Returns the exponential of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the exponential.
 * @returns
 */
function IMEXP(inumber) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Return exponential of complex number
  const e = Math.exp(x);
  return COMPLEX(e * Math.cos(y), e * Math.sin(y), unit)
}

/**
 * Returns the natural logarithm of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the natural logarithm.
 * @returns
 */
function IMLN(inumber) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Return exponential of complex number
  return COMPLEX(Math.log(Math.sqrt(x * x + y * y)), Math.atan(y / x), unit)
}

/**
 * Returns the base-10 logarithm of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the common logarithm.
 * @returns
 */
function IMLOG10(inumber) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Return exponential of complex number
  return COMPLEX(Math.log(Math.sqrt(x * x + y * y)) / Math.log(10), Math.atan(y / x) / Math.log(10), unit)
}

/**
 * Returns the base-2 logarithm of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the base-2 logarithm.
 * @returns
 */
function IMLOG2(inumber) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Return exponential of complex number
  return COMPLEX(Math.log(Math.sqrt(x * x + y * y)) / Math.log(2), Math.atan(y / x) / Math.log(2), unit)
}

/**
 * Returns a complex number raised to an integer power.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number you want to raise to a power.
 * @param {*} number The power to which you want to raise the complex number.
 * @returns
 */
function IMPOWER(inumber, number) {
  number = parseNumber(number);
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(number, x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Calculate power of modulus
  const p = Math.pow(IMABS(inumber), number);

  // Calculate argument
  const t = IMARGUMENT(inumber);

  // Return exponential of complex number
  return COMPLEX(p * Math.cos(number * t), p * Math.sin(number * t), unit)
}

/**
 * Returns the product of complex numbers.
 *
 * Category: Engineering
 *
 * @param {*} args inumber1, [inumber2], … Inumber1 is required, subsequent inumbers are not. 1 to 255 complex numbers to multiply.
 * @returns
 */
function IMPRODUCT() {
  // Initialize result
  let result = arguments[0];

  if (!arguments.length) {
    return value
  }

  // Loop on all numbers
  for (let i = 1; i < arguments.length; i++) {
    // Lookup coefficients of two complex numbers
    const a = IMREAL(result);
    const b = IMAGINARY(result);
    const c = IMREAL(arguments[i]);
    const d = IMAGINARY(arguments[i]);

    if (anyIsError(a, b, c, d)) {
      return value
    }

    // Complute product of two complex numbers
    result = COMPLEX(a * c - b * d, a * d + b * c);
  }

  // Return product of complex numbers
  return result
}

/**
 * Returns the real coefficient of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the real coefficient.
 * @returns
 */
function IMREAL(inumber) {
  if (inumber === undefined || inumber === true || inumber === false) {
    return value
  }

  // Return 0 if inumber is equal to 0
  if (inumber === 0 || inumber === '0') {
    return 0
  }

  // Handle special cases
  if (['i', '+i', '1i', '+1i', '-i', '-1i', 'j', '+j', '1j', '+1j', '-j', '-1j'].indexOf(inumber) >= 0) {
    return 0
  }

  // Force String type
  inumber = inumber + '';

  // Lookup sign
  let plus = inumber.indexOf('+');
  let minus = inumber.indexOf('-');

  if (plus === 0) {
    plus = inumber.indexOf('+', 1);
  }

  if (minus === 0) {
    minus = inumber.indexOf('-', 1);
  }

  // Lookup imaginary unit
  const last = inumber.substring(inumber.length - 1, inumber.length);
  const unit = last === 'i' || last === 'j';

  if (plus >= 0 || minus >= 0) {
    // Return error if imaginary unit is neither i nor j
    if (!unit) {
      return num
    }

    // Return real coefficient of complex number
    if (plus >= 0) {
      return isNaN(inumber.substring(0, plus)) || isNaN(inumber.substring(plus + 1, inumber.length - 1))
        ? num
        : Number(inumber.substring(0, plus))
    } else {
      return isNaN(inumber.substring(0, minus)) || isNaN(inumber.substring(minus + 1, inumber.length - 1))
        ? num
        : Number(inumber.substring(0, minus))
    }
  } else {
    if (unit) {
      return isNaN(inumber.substring(0, inumber.length - 1)) ? num : 0
    } else {
      return isNaN(inumber) ? num : inumber
    }
  }
}

/**
 * Returns the secant of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the secant.
 * @returns
 */
function IMSEC(inumber) {
  // Return error if inumber is a logical value
  if (inumber === true || inumber === false) {
    return value
  }

  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Return secant of complex number
  return IMDIV('1', IMCOS(inumber))
}

/**
 * Returns the hyperbolic secant of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the hyperbolic secant.
 * @returns
 */
function IMSECH(inumber) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Return hyperbolic secant of complex number
  return IMDIV('1', IMCOSH(inumber))
}

/**
 * Returns the sine of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the sine.
 * @returns
 */
function IMSIN(inumber) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Return sine of complex number
  return COMPLEX(
    (Math.sin(x) * (Math.exp(y) + Math.exp(-y))) / 2,
    (Math.cos(x) * (Math.exp(y) - Math.exp(-y))) / 2,
    unit
  )
}

/**
 * Returns the hyperbolic sine of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the hyperbolic sine.
 * @returns
 */
function IMSINH(inumber) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Return hyperbolic sine of complex number
  return COMPLEX(
    (Math.cos(y) * (Math.exp(x) - Math.exp(-x))) / 2,
    (Math.sin(y) * (Math.exp(x) + Math.exp(-x))) / 2,
    unit
  )
}

/**
 * Returns the square root of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the square root.
 * @returns
 */
function IMSQRT(inumber) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Lookup imaginary unit
  let unit = inumber.substring(inumber.length - 1);
  unit = unit === 'i' || unit === 'j' ? unit : 'i';

  // Calculate power of modulus
  const s = Math.sqrt(IMABS(inumber));

  // Calculate argument
  const t = IMARGUMENT(inumber);

  // Return exponential of complex number
  return COMPLEX(s * Math.cos(t / 2), s * Math.sin(t / 2), unit)
}

/**
 * Returns the cosecant of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the cosecant.
 * @returns
 */
function IMCSC(inumber) {
  // Return error if inumber is a logical value
  if (inumber === true || inumber === false) {
    return value
  }

  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  // Return error if either coefficient is not a number
  if (anyIsError(x, y)) {
    return num
  }

  // Return cosecant of complex number
  return IMDIV('1', IMSIN(inumber))
}

/**
 * Returns the hyperbolic cosecant of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the hyperbolic cosecant.
 * @returns
 */
function IMCSCH(inumber) {
  // Return error if inumber is a logical value
  if (inumber === true || inumber === false) {
    return value
  }

  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  // Return error if either coefficient is not a number
  if (anyIsError(x, y)) {
    return num
  }

  // Return hyperbolic cosecant of complex number
  return IMDIV('1', IMSINH(inumber))
}

/**
 * Returns the difference between two complex numbers.
 *
 * Category: Engineering
 *
 * @param {*} inumber1 The complex number from which to subtract inumber2.
 * @param {*} inumber2 The complex number to subtract from inumber1.
 * @returns
 */
function IMSUB(inumber1, inumber2) {
  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const a = IMREAL(inumber1);
  const b = IMAGINARY(inumber1);
  const c = IMREAL(inumber2);
  const d = IMAGINARY(inumber2);

  if (anyIsError(a, b, c, d)) {
    return value
  }

  // Lookup imaginary unit
  const unit1 = inumber1.substring(inumber1.length - 1);
  const unit2 = inumber2.substring(inumber2.length - 1);
  let unit = 'i';

  if (unit1 === 'j') {
    unit = 'j';
  } else if (unit2 === 'j') {
    unit = 'j';
  }

  // Return _ of two complex numbers
  return COMPLEX(a - c, b - d, unit)
}

/**
 * Returns the sum of complex numbers.
 *
 * Category: Engineering
 *
 * @param {*} args inumber1, [inumber2], ... Inumber1 is required, subsequent numbers are not. 1 to 255 complex numbers to add.
 * @returns
 */
function IMSUM() {
  if (!arguments.length) {
    return value
  }

  const args = flatten(arguments);

  // Initialize result
  let result = args[0];

  // Loop on all numbers
  for (let i = 1; i < args.length; i++) {
    // Lookup coefficients of two complex numbers
    const a = IMREAL(result);
    const b = IMAGINARY(result);
    const c = IMREAL(args[i]);
    const d = IMAGINARY(args[i]);

    if (anyIsError(a, b, c, d)) {
      return value
    }

    // Complute product of two complex numbers
    result = COMPLEX(a + c, b + d);
  }

  // Return sum of complex numbers
  return result
}

/**
 * Returns the tangent of a complex number.
 *
 * Category: Engineering
 *
 * @param {*} inumber A complex number for which you want the tangent.
 * @returns
 */
function IMTAN(inumber) {
  // Return error if inumber is a logical value
  if (inumber === true || inumber === false) {
    return value
  }

  // Lookup real and imaginary coefficients using Formula.js [http://formulajs.org]
  const x = IMREAL(inumber);
  const y = IMAGINARY(inumber);

  if (anyIsError(x, y)) {
    return value
  }

  // Return tangent of complex number
  return IMDIV(IMSIN(inumber), IMCOS(inumber))
}

/**
 * Converts an octal number to binary.
 *
 * Category: Engineering
 *
 * @param {*} number The octal number you want to convert. Number may not contain more than 10 characters. The most significant bit of number is the sign bit. The remaining 29 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @param {*} places Optional. The number of characters to use. If places is omitted, OCT2BIN uses the minimum number of characters necessary. Places is useful for padding the return value with leading 0s (zeros).
 * @returns
 */
function OCT2BIN(number, places) {
  // Return error if number is not hexadecimal or contains more than ten characters (10 digits)
  if (!/^[0-7]{1,10}$/.test(number)) {
    return num
  }

  // Check if number is negative
  const negative = !!(number.length === 10 && number.substring(0, 1) === '7');

  // Convert octal number to decimal
  const decimal = negative ? parseInt(number, 8) - 1073741824 : parseInt(number, 8);

  // Return error if number is lower than -512 or greater than 511
  if (decimal < -512 || decimal > 511) {
    return num
  }

  // Ignore places and return a 10-character binary number if number is negative
  if (negative) {
    return '1' + REPT('0', 9 - (512 + decimal).toString(2).length) + (512 + decimal).toString(2)
  }

  // Convert decimal number to binary
  const result = decimal.toString(2);

  // Return binary number using the minimum number of characters necessary if places is undefined
  if (typeof places === 'undefined') {
    return result
  } else {
    // Return error if places is nonnumeric
    if (isNaN(places)) {
      return value
    }

    // Return error if places is negative
    if (places < 0) {
      return num
    }

    // Truncate places in case it is not an integer
    places = Math.floor(places);

    // Pad return value with leading 0s (zeros) if necessary (using Underscore.string)
    return places >= result.length ? REPT('0', places - result.length) + result : num
  }
}

/**
 * Converts an octal number to decimal.
 *
 * Category: Engineering
 *
 * @param {*} number The octal number you want to convert. Number may not contain more than 10 octal characters (30 bits). The most significant bit of number is the sign bit. The remaining 29 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @returns
 */
function OCT2DEC(number) {
  // Return error if number is not octal or contains more than ten characters (10 digits)
  if (!/^[0-7]{1,10}$/.test(number)) {
    return num
  }

  // Convert octal number to decimal
  const decimal = parseInt(number, 8);

  // Return decimal number
  return decimal >= 536870912 ? decimal - 1073741824 : decimal
}

/**
 * Converts an octal number to hexadecimal.
 *
 * Category: Engineering
 *
 * @param {*} number The octal number you want to convert. Number may not contain more than 10 octal characters (30 bits). The most significant bit of number is the sign bit. The remaining 29 bits are magnitude bits. Negative numbers are represented using two's-complement notation.
 * @param {*} places Optional. The number of characters to use. If places is omitted, OCT2HEX uses the minimum number of characters necessary. Places is useful for padding the return value with leading 0s (zeros).
 * @returns
 */
function OCT2HEX(number, places) {
  // Return error if number is not octal or contains more than ten characters (10 digits)
  if (!/^[0-7]{1,10}$/.test(number)) {
    return num
  }

  // Convert octal number to decimal
  const decimal = parseInt(number, 8);

  // Ignore places and return a 10-character octal number if number is negative
  if (decimal >= 536870912) {
    return 'ff' + (decimal + 3221225472).toString(16)
  }

  // Convert decimal number to hexadecimal
  const result = decimal.toString(16);

  // Return hexadecimal number using the minimum number of characters necessary if places is undefined
  if (places === undefined) {
    return result
  } else {
    // Return error if places is nonnumeric
    if (isNaN(places)) {
      return value
    }

    // Return error if places is negative
    if (places < 0) {
      return num
    }

    // Truncate places in case it is not an integer
    places = Math.floor(places);

    // Pad return value with leading 0s (zeros) if necessary (using Underscore.string)
    return places >= result.length ? REPT('0', places - result.length) + result : num
  }
}

const BETADIST = BETA.DIST;
const BETAINV = BETA.INV;
const BINOMDIST = BINOM.DIST;
const CEILINGMATH = CEILING.MATH;
const CEILINGPRECISE = CEILING.PRECISE;
const CHIDIST = CHISQ.DIST;
const CHIDISTRT = CHISQ.DIST.RT;
const CHIINV = CHISQ.INV;
const CHIINVRT = CHISQ.INV.RT;
const CHITEST = CHISQ.TEST;
const COVAR = COVARIANCE.P;
const COVARIANCEP = COVARIANCE.P;
const COVARIANCES = COVARIANCE.S;
const CRITBINOM = BINOM.INV;
const ERFCPRECISE = ERFC.PRECISE;
const ERFPRECISE = ERF.PRECISE;
const EXPONDIST = EXPON.DIST;
const FDIST = F.DIST;
const FDISTRT = F.DIST.RT;
const FINV = F.INV;
const FINVRT = F.INV.RT;
const FLOORMATH = FLOOR.MATH;
const FLOORPRECISE = FLOOR.PRECISE;
const FTEST = F.TEST;
const GAMMADIST = GAMMA.DIST;
const GAMMAINV = GAMMA.INV;
const GAMMALNPRECISE = GAMMALN.PRECISE;
const HYPGEOMDIST = HYPGEOM.DIST;
const LOGINV = LOGNORM.INV;
const LOGNORMDIST = LOGNORM.DIST;
const LOGNORMINV = LOGNORM.INV;
const MODEMULT = MODE.MULT;
const MODESNGL = MODE.SNGL;
const NEGBINOMDIST = NEGBINOM.DIST;
const NETWORKDAYSINTL = NETWORKDAYS.INTL;
const NORMDIST = NORM.DIST;
const NORMINV = NORM.INV;
const NORMSDIST = NORM.S.DIST;
const NORMSINV = NORM.S.INV;
const PERCENTILEEXC = PERCENTILE.EXC;
const PERCENTILEINC = PERCENTILE.INC;
const PERCENTRANKEXC = PERCENTRANK.EXC;
const PERCENTRANKINC = PERCENTRANK.INC;
const POISSONDIST = POISSON.DIST;
const QUARTILEEXC = QUARTILE.EXC;
const QUARTILEINC = QUARTILE.INC;
const RANKAVG = RANK.AVG;
const RANKEQ = RANK.EQ;
const SKEWP = SKEW.P;
const STDEVP = STDEV.P;
const STDEVS = STDEV.S;
const TDIST = T$1.DIST;
const TDISTRT = T$1.DIST.RT;
const TINV = T$1.INV;
const TTEST = T$1.TEST;
const VARP = VAR.P;
const VARS = VAR.S;
const WEIBULLDIST = WEIBULL.DIST;
const WORKDAYINTL = WORKDAY.INTL;
const ZTEST = Z.TEST;

function compact(array) {
  const result = [];

  arrayEach(array, (value) => {
    if (value) {
      result.push(value);
    }
  });

  return result
}

function findResultIndex(database, criterias) {
  const matches = {};

  for (let i = 1; i < database[0].length; ++i) {
    matches[i] = true;
  }

  let maxCriteriaLength = criterias[0].length;

  for (let i = 1; i < criterias.length; ++i) {
    if (criterias[i].length > maxCriteriaLength) {
      maxCriteriaLength = criterias[i].length;
    }
  }

  for (let k = 1; k < database.length; ++k) {
    for (let l = 1; l < database[k].length; ++l) {
      let currentCriteriaResult = false;
      let hasMatchingCriteria = false;

      for (let j = 0; j < criterias.length; ++j) {
        const criteria = criterias[j];

        if (criteria.length < maxCriteriaLength) {
          continue
        }

        const criteriaField = criteria[0];

        if (database[k][0] !== criteriaField) {
          continue
        }

        hasMatchingCriteria = true;

        for (let p = 1; p < criteria.length; ++p) {
          if (!currentCriteriaResult) {
            const isWildcard = criteria[p] === void 0 || criteria[p] === '*';

            if (isWildcard) {
              currentCriteriaResult = true;
            } else {
              const tokenizedCriteria = parse(criteria[p] + '');
              const tokens = [createToken(database[k][l], TOKEN_TYPE_LITERAL)].concat(
                tokenizedCriteria
              );

              currentCriteriaResult = compute(tokens);
            }
          }
        }
      }

      if (hasMatchingCriteria) {
        matches[l] = matches[l] && currentCriteriaResult;
      }
    }
  }

  const result = [];

  for (let n = 0; n < database[0].length; ++n) {
    if (matches[n]) {
      result.push(n - 1);
    }
  }

  return result
}

// Database functions
/**
 * Returns the average of selected database entries.
 *
 * Category: Database
 *
 * @param {*} database Range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria Range of values that contains the conditions you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DAVERAGE(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  let sum = 0;

  arrayEach(resultIndexes, (value) => {
    sum += targetFields[value];
  });

  return resultIndexes.length === 0 ? div0 : sum / resultIndexes.length
}

/**
 * Counts the values that contain numbers in a database.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as the argument includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DCOUNT(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  const targetValues = [];

  arrayEach(resultIndexes, (value) => {
    targetValues.push(targetFields[value]);
  });

  return COUNT(targetValues)
}

/**
 * Counts nonblank values in a database.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Optional. Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DCOUNTA(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  const targetValues = [];

  arrayEach(resultIndexes, (value) => {
    targetValues.push(targetFields[value]);
  });

  return COUNTA(targetValues)
}

/**
 * Extracts from a database a single record that matches the specified criteria.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DGET(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  // Return error if no record meets the criteria
  if (resultIndexes.length === 0) {
    return value
  }
  // Returns the #NUM! error value because more than one record meets the
  // criteria
  if (resultIndexes.length > 1) {
    return num
  }

  return targetFields[resultIndexes[0]]
}

/**
 * Returns the maximum value from selected database entries.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DMAX(database, field, criteria) {
  // Return error if field is not a number and not a string

  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  let maxValue = targetFields[resultIndexes[0]];

  arrayEach(resultIndexes, (value) => {
    if (maxValue < targetFields[value]) {
      maxValue = targetFields[value];
    }
  });

  return maxValue
}

/**
 * Returns the minimum value from selected database entries.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DMIN(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  let minValue = targetFields[resultIndexes[0]];

  arrayEach(resultIndexes, (value) => {
    if (minValue > targetFields[value]) {
      minValue = targetFields[value];
    }
  });

  return minValue
}

/**
 * Multiplies the values in a particular field of records that match the criteria in a database.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DPRODUCT(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  let targetValues = [];

  arrayEach(resultIndexes, (value) => {
    targetValues.push(targetFields[value]);
  });
  targetValues = compact(targetValues);

  let result = 1;

  arrayEach(targetValues, (value) => {
    result *= value;
  });

  return result
}

/**
 * Estimates the standard deviation based on a sample of selected database entries.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DSTDEV(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  let targetValues = [];

  arrayEach(resultIndexes, (value) => {
    targetValues.push(targetFields[value]);
  });

  targetValues = compact(targetValues);

  return STDEV.S(targetValues)
}

/**
 * Calculates the standard deviation based on the entire population of selected database entries.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DSTDEVP(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  let targetValues = [];

  arrayEach(resultIndexes, (value) => {
    targetValues.push(targetFields[value]);
  });

  targetValues = compact(targetValues);

  return STDEV.P(targetValues)
}

/**
 * Adds the numbers in the field column of records in the database that match the criteria.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria Is the range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DSUM(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  const targetValues = [];

  arrayEach(resultIndexes, (value) => {
    targetValues.push(targetFields[value]);
  });

  return SUM(targetValues)
}

/**
 * Estimates variance based on a sample from selected database entries.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DVAR(database, field, criteria) {
  // Return error if field is not a number and not a string
  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  const targetValues = [];

  arrayEach(resultIndexes, (value) => {
    targetValues.push(targetFields[value]);
  });

  return VAR.S(targetValues)
}

/**
 * Calculates variance based on the entire population of selected database entries.
 *
 * Category: Database
 *
 * @param {*} database The range of values that makes up the list or database. A database is a list of related data in which rows of related information are records, and columns of data are fields. The first row of the list contains labels for each column.
 * @param {*} field Indicates which column is used in the function. Enter the column label enclosed between double quotation marks, such as "Age" or "Yield," or a number (without quotation marks) that represents the position of the column within the list: 1 for the first column, 2 for the second column, and so on.
 * @param {*} criteria The range of values that contains the conditions that you specify. You can use any range for the criteria argument, as long as it includes at least one column label and at least one value below the column label in which you specify a condition for the column.
 * @returns
 */
function DVARP(database, field, criteria) {
  // Return error if field is not a number and not a string

  if (isNaN(field) && typeof field !== 'string') {
    return value
  }

  const resultIndexes = findResultIndex(database, criteria);
  let targetFields = [];

  if (typeof field === 'string') {
    const index = findField(database, field);
    targetFields = rest(database[index]);
  } else {
    targetFields = rest(database[field]);
  }

  const targetValues = [];

  arrayEach(resultIndexes, (value) => {
    targetValues.push(targetFields[value]);
  });

  return VAR.P(targetValues)
}

function validDate(d) {
  return d && d.getTime && !isNaN(d.getTime())
}

function ensureDate(d) {
  return d instanceof Date ? d : new Date(d)
}

/**
 * Returns the accrued interest for a security that pays periodic interest.
 *
 * Category: Financial
 *
 * @param {*} issue The security's issue date.
 * @param {*} first_interest The security's first interest date.
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} rate The security's annual coupon rate.
 * @param {*} par The security's par value. If you omit par, ACCRINT uses $1,000.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @param {*} calc_method Optional. Not implemented in formulajs. A logical value that specifies the way to calculate the total accrued interest when the date of settlement is later than the date of first_interest. A value of TRUE (1) returns the total accrued interest from issue to settlement. A value of FALSE (0) returns the accrued interest from first_interest to settlement. If you do not enter the argument, it defaults to TRUE.
 * @returns
 */
function ACCRINT(issue, first_interest, settlement, rate, par, frequency, basis) {
  // Return error if either date is invalid
  issue = ensureDate(issue);
  first_interest = ensureDate(first_interest);
  settlement = ensureDate(settlement);

  if (!validDate(issue) || !validDate(first_interest) || !validDate(settlement)) {
    return value
  }

  // Return error if either rate or par are lower than or equal to zero
  if (rate <= 0 || par <= 0) {
    return num
  }

  // Return error if frequency is neither 1, 2, or 4
  if ([1, 2, 4].indexOf(frequency) === -1) {
    return num
  }

  // Return error if basis is neither 0, 1, 2, 3, or 4
  if ([0, 1, 2, 3, 4].indexOf(basis) === -1) {
    return num
  }

  // Return error if settlement is before or equal to issue
  if (settlement <= issue) {
    return num
  }

  // Set default values
  par = par || 0;
  basis = basis || 0;

  // Compute accrued interest
  return par * rate * YEARFRAC(issue, settlement, basis)
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the accrued interest for a security that pays interest at maturity.
 *
 * Category: Financial
 *
 * @param {*} issue The security's issue date.
 * @param {*} settlement The security's maturity date.
 * @param {*} rate The security's annual coupon rate.
 * @param {*} par The security's par value. If you omit par, ACCRINTM uses $1,000.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function ACCRINTM() {
  throw new Error('ACCRINTM is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the depreciation for each accounting period by using a depreciation coefficient.
 *
 * Category: Financial
 *
 * @param {*} cost The cost of the asset.
 * @param {*} date_purchased The date of the purchase of the asset.
 * @param {*} first_period The date of the end of the first period.
 * @param {*} salvage The salvage value at the end of the life of the asset.
 * @param {*} period The period.
 * @param {*} rate The rate of depreciation.
 * @param {*} basis Optional. The year basis to be used.
 * @returns
 */
function AMORDEGRC() {
  throw new Error('AMORDEGRC is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the depreciation for each accounting period.
 *
 * Category: Financial
 *
 * @param {*} cost The cost of the asset.
 * @param {*} date_purchased The date of the purchase of the asset.
 * @param {*} first_period The date of the end of the first period.
 * @param {*} salvage The salvage value at the end of the life of the asset.
 * @param {*} period The period.
 * @param {*} rate The rate of depreciation.
 * @param {*} basis Optional. The year basis to be used.
 * @returns
 */
function AMORLINC() {
  throw new Error('AMORLINC is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the number of days from the beginning of the coupon period to the settlement date.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function COUPDAYBS() {
  throw new Error('COUPDAYBS is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the number of days in the coupon period that contains the settlement date.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function COUPDAYS() {
  throw new Error('COUPDAYS is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the number of days from the settlement date to the next coupon date.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function COUPDAYSNC() {
  throw new Error('COUPDAYSNC is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the next coupon date after the settlement date.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function COUPNCD() {
  throw new Error('COUPNCD is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the number of coupons payable between the settlement date and maturity date.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function COUPNUM() {
  throw new Error('COUPNUM is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the previous coupon date before the settlement date.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function COUPPCD() {
  throw new Error('COUPPCD is not implemented')
}

/**
 * Returns the cumulative interest paid between two periods.
 *
 * Category: Financial
 *
 * @param {*} rate The interest rate.
 * @param {*} nper The total number of payment periods.
 * @param {*} pv The present value.
 * @param {*} start_period The first period in the calculation. Payment periods are numbered beginning with 1.
 * @param {*} end_period The last period in the calculation.
 * @param {*} type The timing of the payment.
 * @returns
 */
function CUMIPMT(rate, nper, pv, start_period, end_period, type) {
  rate = parseNumber(rate);
  nper = parseNumber(nper);
  pv = parseNumber(pv);

  if (anyIsError(rate, nper, pv)) {
    return value
  }

  if (rate <= 0 || nper <= 0 || pv <= 0) {
    return num
  }

  if (start_period < 1 || end_period < 1 || start_period > end_period) {
    return num
  }

  if (type !== 0 && type !== 1) {
    return num
  }

  const payment = PMT(rate, nper, pv, 0, type);
  let interest = 0;

  if (start_period === 1) {
    if (type === 0) {
      interest = -pv;
    }

    start_period++;
  }

  for (let i = start_period; i <= end_period; i++) {
    interest += type === 1 ? FV(rate, i - 2, payment, pv, 1) - payment : FV(rate, i - 1, payment, pv, 0);
  }

  interest *= rate;

  return interest
}

/**
 * Returns the cumulative principal paid on a loan between two periods.
 *
 * Category: Financial
 *
 * @param {*} rate The interest rate.
 * @param {*} nper The total number of payment periods.
 * @param {*} pv The present value.
 * @param {*} start_period The first period in the calculation. Payment periods are numbered beginning with 1.
 * @param {*} end_period The last period in the calculation.
 * @param {*} type The timing of the payment.
 * @returns
 */
function CUMPRINC(rate, nper, pv, start_period, end, type) {
  // Credits: algorithm inspired by Apache OpenOffice
  // Credits: Hannes Stiebitzhofer for the translations of function and variable names
  rate = parseNumber(rate);
  nper = parseNumber(nper);
  pv = parseNumber(pv);

  if (anyIsError(rate, nper, pv)) {
    return value
  }

  // Return error if either rate, nper, or value are lower than or equal to zero
  if (rate <= 0 || nper <= 0 || pv <= 0) {
    return num
  }

  // Return error if start < 1, end < 1, or start > end
  if (start_period < 1 || end < 1 || start_period > end) {
    return num
  }

  // Return error if type is neither 0 nor 1
  if (type !== 0 && type !== 1) {
    return num
  }

  // Compute cumulative principal
  const payment = PMT(rate, nper, pv, 0, type);
  let principal = 0;

  if (start_period === 1) {
    principal = type === 0 ? payment + pv * rate : payment;

    start_period++;
  }

  for (let i = start_period; i <= end; i++) {
    principal +=
      type > 0
        ? payment - (FV(rate, i - 2, payment, pv, 1) - payment) * rate
        : payment - FV(rate, i - 1, payment, pv, 0) * rate;
  }

  // Return cumulative principal
  return principal
}

/**
 * Returns the depreciation of an asset for a specified period by using the fixed-declining balance method.
 *
 * Category: Financial
 *
 * @param {*} cost The initial cost of the asset.
 * @param {*} salvage The value at the end of the depreciation (sometimes called the salvage value of the asset).
 * @param {*} life The number of periods over which the asset is being depreciated (sometimes called the useful life of the asset).
 * @param {*} period The period for which you want to calculate the depreciation. Period must use the same units as life.
 * @param {*} month Optional. The number of months in the first year. If month is omitted, it is assumed to be 12.
 * @returns
 */
function DB(cost, salvage, life, period, month) {
  // Initialize month
  month = month === undefined ? 12 : month;

  cost = parseNumber(cost);
  salvage = parseNumber(salvage);
  life = parseNumber(life);
  period = parseNumber(period);
  month = parseNumber(month);

  if (anyIsError(cost, salvage, life, period, month)) {
    return value
  }

  // Return error if any of the parameters is negative
  if (cost < 0 || salvage < 0 || life < 0 || period < 0) {
    return num
  }

  // Return error if month is not an integer between 1 and 12
  if ([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12].indexOf(month) === -1) {
    return num
  }

  // Return error if period is greater than life
  if (period > life) {
    return num
  }

  // Return 0 (zero) if salvage is greater than or equal to cost
  if (salvage >= cost) {
    return 0
  }

  // Rate is rounded to three decimals places
  const rate = (1 - Math.pow(salvage / cost, 1 / life)).toFixed(3);

  // Compute initial depreciation
  const initial = (cost * rate * month) / 12;

  // Compute total depreciation
  let total = initial;
  let current = 0;
  const ceiling = period === life ? life - 1 : period;

  for (let i = 2; i <= ceiling; i++) {
    current = (cost - total) * rate;
    total += current;
  }

  // Depreciation for the first and last periods are special cases
  if (period === 1) {
    // First period
    return initial
  } else if (period === life) {
    // Last period

    return (cost - total) * rate
  } else {
    return current
  }
}

/**
 * Returns the depreciation of an asset for a specified period by using the double-declining balance method or some other method that you specify.
 *
 * Category: Financial
 *
 * @param {*} cost The initial cost of the asset.
 * @param {*} salvage The value at the end of the depreciation (sometimes called the salvage value of the asset). This value can be 0.
 * @param {*} life The number of periods over which the asset is being depreciated (sometimes called the useful life of the asset).
 * @param {*} period The period for which you want to calculate the depreciation. Period must use the same units as life.
 * @param {*} factor Optional. The rate at which the balance declines. If factor is omitted, it is assumed to be 2 (the double-declining balance method).
 * @returns
 */
function DDB(cost, salvage, life, period, factor) {
  // Initialize factor
  factor = factor === undefined ? 2 : factor;

  cost = parseNumber(cost);
  salvage = parseNumber(salvage);
  life = parseNumber(life);
  period = parseNumber(period);
  factor = parseNumber(factor);

  if (anyIsError(cost, salvage, life, period, factor)) {
    return value
  }

  // Return error if any of the parameters is negative or if factor is null
  if (cost < 0 || salvage < 0 || life < 0 || period < 0 || factor <= 0) {
    return num
  }

  // Return error if period is greater than life
  if (period > life) {
    return num
  }

  // Return 0 (zero) if salvage is greater than or equal to cost
  if (salvage >= cost) {
    return 0
  }

  // Compute depreciation
  let total = 0;
  let current = 0;

  for (let i = 1; i <= period; i++) {
    current = Math.min((cost - total) * (factor / life), cost - salvage - total);
    total += current;
  }

  // Return depreciation
  return current
}

/**
 * Returns the discount rate for a security.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} pr The security's price per $100 face value.
 * @param {*} redemption The security's redemption value per $100 face value.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function DISC(settlement, maturity, pr, redemption, basis) {
  settlement = parseDate(settlement);
  maturity = parseDate(maturity);
  pr = parseNumber(pr);
  redemption = parseNumber(redemption);
  basis = parseNumber(basis);

  basis = basis || 0;

  if (anyIsError(settlement, maturity, pr, redemption, basis)) {
    return value
  }

  if (pr <= 0 || redemption <= 0) {
    return num
  }

  if (settlement >= maturity) {
    return value
  }

  let basisVal, diff;
  switch (basis) {
    case 0:
      basisVal = 360;
      diff = DAYS360(settlement, maturity, false);
      break
    case 1:
      basisVal = 365;
      diff = DATEDIF(settlement, maturity, 'D');
      break
    case 2:
      basisVal = 360;
      diff = DATEDIF(settlement, maturity, 'D');
      break
    case 3:
      basisVal = 365;
      diff = DATEDIF(settlement, maturity, 'D');
      break
    case 4:
      basisVal = 360;
      diff = DAYS360(settlement, maturity, true);
      break
    default:
      return num
  }

  return (((redemption - pr) / redemption) * basisVal) / diff
}

/**
 * Converts a dollar price, expressed as a fraction, into a dollar price, expressed as a decimal number.
 *
 * Category: Financial
 *
 * @param {*} fractional_dollar A number expressed as an integer part and a fraction part, separated by a decimal symbol.
 * @param {*} fraction The integer to use in the denominator of the fraction.
 * @returns
 */
function DOLLARDE(fractional_dollar, fraction) {
  // Credits: algorithm inspired by Apache OpenOffice
  fractional_dollar = parseNumber(fractional_dollar);
  fraction = parseNumber(fraction);

  if (anyIsError(fractional_dollar, fraction)) {
    return value
  }

  // Return error if fraction is negative
  if (fraction < 0) {
    return num
  }

  // Return error if fraction is greater than or equal to 0 and less than 1
  if (fraction >= 0 && fraction < 1) {
    return div0
  }

  // Truncate fraction if it is not an integer
  fraction = parseInt(fraction, 10);

  // Compute integer part
  let result = parseInt(fractional_dollar, 10);

  // Add decimal part
  result += ((fractional_dollar % 1) * Math.pow(10, Math.ceil(Math.log(fraction) / Math.LN10))) / fraction;

  // Round result
  const power = Math.pow(10, Math.ceil(Math.log(fraction) / Math.LN2) + 1);
  result = Math.round(result * power) / power;

  // Return converted dollar price
  return result
}

/**
 * Converts a dollar price, expressed as a decimal number, into a dollar price, expressed as a fraction.
 *
 * Category: Financial
 *
 * @param {*} decimal_dollar A decimal number.
 * @param {*} fraction The integer to use in the denominator of a fraction.
 * @returns
 */
function DOLLARFR(decimal_dollar, fraction) {
  // Credits: algorithm inspired by Apache OpenOffice
  decimal_dollar = parseNumber(decimal_dollar);
  fraction = parseNumber(fraction);

  if (anyIsError(decimal_dollar, fraction)) {
    return value
  }

  // Return error if fraction is negative
  if (fraction < 0) {
    return num
  }

  // Return error if fraction is greater than or equal to 0 and less than 1
  if (fraction >= 0 && fraction < 1) {
    return div0
  }

  // Truncate fraction if it is not an integer
  fraction = parseInt(fraction, 10);

  // Compute integer part
  let result = parseInt(decimal_dollar, 10);

  // Add decimal part
  result += (decimal_dollar % 1) * Math.pow(10, -Math.ceil(Math.log(fraction) / Math.LN10)) * fraction;

  // Return converted dollar price
  return result
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the annual duration of a security with periodic interest payments.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} coupon The security's annual coupon rate.
 * @param {*} yld The security's annual yield.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function DURATION() {
  throw new Error('DURATION is not implemented')
}

/**
 * Returns the effective annual interest rate.
 *
 * Category: Financial
 *
 * @param {*} nominal_rate The nominal interest rate.
 * @param {*} npery The number of compounding periods per year.
 * @returns
 */
function EFFECT(nominal_rate, npery) {
  nominal_rate = parseNumber(nominal_rate);
  npery = parseNumber(npery);

  if (anyIsError(nominal_rate, npery)) {
    return value
  }

  // Return error if rate <=0 or periods < 1
  if (nominal_rate <= 0 || npery < 1) {
    return num
  }

  // Truncate periods if it is not an integer
  npery = parseInt(npery, 10);

  // Return effective annual interest rate
  return Math.pow(1 + nominal_rate / npery, npery) - 1
}

/**
 * Returns the future value of an investment.
 *
 * Category: Financial
 *
 * @param {*} rate The interest rate per period.
 * @param {*} nper The total number of payment periods in an annuity.
 * @param {*} pmt The payment made each period; it cannot change over the life of the annuity. Typically, pmt contains principal and interest but no other fees or taxes. If pmt is omitted, you must include the pv argument.
 * @param {*} pv Optional. The present value, or the lump-sum amount that a series of future payments is worth right now. If pv is omitted, it is assumed to be 0 (zero), and you must include the pmt argument.
 * @param {*} type Optional. The number 0 or 1 and indicates when payments are due. If type is omitted, it is assumed to be 0.
 * @returns
 */
function FV(rate, nper, payment, value$1, type) {
  // Credits: algorithm inspired by Apache OpenOffice
  value$1 = value$1 || 0;
  type = type || 0;

  rate = parseNumber(rate);
  nper = parseNumber(nper);
  payment = parseNumber(payment);
  value$1 = parseNumber(value$1);
  type = parseNumber(type);

  if (anyIsError(rate, nper, payment, value$1, type)) {
    return value
  }

  // Return future value
  let result;

  if (rate === 0) {
    result = value$1 + payment * nper;
  } else {
    const term = Math.pow(1 + rate, nper);

    result =
      type === 1
        ? value$1 * term + (payment * (1 + rate) * (term - 1)) / rate
        : value$1 * term + (payment * (term - 1)) / rate;
  }

  return -result
}

/**
 * Returns the future value of an initial principal after applying a series of compound interest rates.
 *
 * Category: Financial
 *
 * @param {*} principal The present value.
 * @param {*} schedule An array of interest rates to apply.
 * @returns
 */
function FVSCHEDULE(principal, schedule) {
  principal = parseNumber(principal);
  schedule = parseNumberArray(flatten(schedule));

  if (anyIsError(principal, schedule)) {
    return value
  }

  const n = schedule.length;
  let future = principal;

  // Apply all interests in schedule

  for (let i = 0; i < n; i++) {
    // Apply scheduled interest
    future *= 1 + schedule[i];
  }

  // Return future value
  return future
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the interest rate for a fully invested security.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} investment The amount invested in the security.
 * @param {*} redemption The amount to be received at maturity.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function INTRATE() {
  throw new Error('INTRATE is not implemented')
}

/**
 * Returns the interest payment for an investment for a given period.
 *
 * Category: Financial
 *
 * @param {*} rate The interest rate per period.
 * @param {*} per The period for which you want to find the interest and must be in the range 1 to nper.
 * @param {*} nper The total number of payment periods in an annuity.
 * @param {*} pv The present value, or the lump-sum amount that a series of future payments is worth right now.
 * @param {*} fv Optional. The future value, or a cash balance you want to attain after the last payment is made. If fv is omitted, it is assumed to be 0 (the future value of a loan, for example, is 0).
 * @param {*} type Optional. The number 0 or 1 and indicates when payments are due. If type is omitted, it is assumed to be 0.
 * @returns
 */
function IPMT(rate, per, nper, pv, fv, type) {
  // Credits: algorithm inspired by Apache OpenOffice
  fv = fv || 0;
  type = type || 0;

  rate = parseNumber(rate);
  per = parseNumber(per);
  nper = parseNumber(nper);
  pv = parseNumber(pv);
  fv = parseNumber(fv);
  type = parseNumber(type);

  if (anyIsError(rate, per, nper, pv, fv, type)) {
    return value
  }

  // Compute payment
  const payment = PMT(rate, nper, pv, fv, type);

  // Compute interest
  let interest =
    per === 1
      ? type === 1
        ? 0
        : -pv
      : type === 1
      ? FV(rate, per - 2, payment, pv, 1) - payment
      : FV(rate, per - 1, payment, pv, 0);

  // Return interest
  return interest * rate
}

/**
 * Returns the internal rate of return for a series of cash flows.
 *
 * Category: Financial
 *
 * @param {*} values An array or a reference to values that contain numbers for which you want to calculate the internal rate of return.
 - Values must contain at least one positive value and one negative value to calculate the internal rate of return.
 - IRR uses the order of values to interpret the order of cash flows. Be sure to enter your payment and income values in the sequence you want.
 - If an array or reference argument contains text, logical values, or empty values, those values are ignored.
 * @param {*} guess Optional. A number that you guess is close to the result of IRR.
 - Microsoft Excel uses an iterative technique for calculating IRR. Starting with guess, IRR cycles through the calculation until the result is accurate within 0.00001 percent. If IRR can't find a result that works after 20 tries, the #NUM! error value is returned.
 - In most cases you do not need to provide guess for the IRR calculation. If guess is omitted, it is assumed to be 0.1 (10 percent).
 - If IRR gives the #NUM! error value, or if the result is not close to what you expected, try again with a different value for guess.
 * @returns
 */
function IRR(values, guess) {
  // Credits: algorithm inspired by Apache OpenOffice
  guess = guess || 0;

  values = parseNumberArray(flatten(values));
  guess = parseNumber(guess);

  if (anyIsError(values, guess)) {
    return value
  }

  // Calculates the resulting amount
  const irrResult = (values, dates, rate) => {
    const r = rate + 1;
    let result = values[0];

    for (let i = 1; i < values.length; i++) {
      result += values[i] / Math.pow(r, (dates[i] - dates[0]) / 365);
    }

    return result
  };

  // Calculates the first derivation
  const irrResultDeriv = (values, dates, rate) => {
    const r = rate + 1;
    let result = 0;

    for (let i = 1; i < values.length; i++) {
      const frac = (dates[i] - dates[0]) / 365;
      result -= (frac * values[i]) / Math.pow(r, frac + 1);
    }

    return result
  };

  // Initialize dates and check that values contains at least one positive value and one negative value
  const dates = [];
  let positive = false;
  let negative = false;

  for (let i = 0; i < values.length; i++) {
    dates[i] = i === 0 ? 0 : dates[i - 1] + 365;

    if (values[i] > 0) {
      positive = true;
    }

    if (values[i] < 0) {
      negative = true;
    }
  }

  // Return error if values does not contain at least one positive value and one negative value
  if (!positive || !negative) {
    return num
  }

  // Initialize guess and resultRate
  guess = guess === undefined ? 0.1 : guess;
  let resultRate = guess;

  // Set maximum epsilon for end of iteration
  const epsMax = 1e-10;

  // Implement Newton's method
  let newRate, epsRate, resultValue;
  let contLoop = true;
  do {
    resultValue = irrResult(values, dates, resultRate);
    newRate = resultRate - resultValue / irrResultDeriv(values, dates, resultRate);
    epsRate = Math.abs(newRate - resultRate);
    resultRate = newRate;
    contLoop = epsRate > epsMax && Math.abs(resultValue) > epsMax;
  } while (contLoop)

  // Return internal rate of return
  return resultRate
}

/**
 * Calculates the interest paid during a specific period of an investment.
 *
 * Category: Financial
 *
 * @param {*} rate The interest rate for the investment.
 * @param {*} per The period for which you want to find the interest, and must be between 1 and Nper.
 * @param {*} nper The total number of payment periods for the investment.
 * @param {*} pv The present value of the investment. For a loan, Pv is the loan amount.
 *
 * @returns
 */
function ISPMT(rate, per, nper, pv) {
  rate = parseNumber(rate);
  per = parseNumber(per);
  nper = parseNumber(nper);
  pv = parseNumber(pv);

  if (anyIsError(rate, per, nper, pv)) {
    return value
  }

  // Return interest
  return pv * rate * (per / nper - 1)
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the Macauley modified duration for a security with an assumed par value of $100.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} coupon The security's annual coupon rate.
 * @param {*} yld The security's annual yield.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function MDURATION() {
  throw new Error('MDURATION is not implemented')
}

/**
 * Returns the internal rate of return where positive and negative cash flows are financed at different rates.
 *
 * Category: Financial
 *
 * @param {*} values An array or a reference to values that contain numbers. These numbers represent a series of payments (negative values) and income (positive values) occurring at regular periods.
 - Values must contain at least one positive value and one negative value to calculate the modified internal rate of return. Otherwise, MIRR returns the #DIV/0! error value.
 - If an array or reference argument contains text, logical values, or empty values, those values are ignored; however, values with the value zero are included.
 * @param {*} finance_rate The interest rate you pay on the money used in the cash flows.
 * @param {*} reinvest_rate The interest rate you receive on the cash flows as you reinvest them.
 * @returns
 */
function MIRR(values, finance_rate, reinvest_rate) {
  values = parseNumberArray(flatten(values));
  finance_rate = parseNumber(finance_rate);
  reinvest_rate = parseNumber(reinvest_rate);

  if (anyIsError(values, finance_rate, reinvest_rate)) {
    return value
  }

  // Initialize number of values
  const n = values.length;

  // Lookup payments (negative values) and incomes (positive values)
  const payments = [];
  const incomes = [];

  for (let i = 0; i < n; i++) {
    if (values[i] < 0) {
      payments.push(values[i]);
    } else {
      incomes.push(values[i]);
    }
  }

  // Return modified internal rate of return
  const num = -NPV(reinvest_rate, incomes) * Math.pow(1 + reinvest_rate, n - 1);
  const den = NPV(finance_rate, payments) * (1 + finance_rate);

  return Math.pow(num / den, 1 / (n - 1)) - 1
}

/**
 * Returns the annual nominal interest rate.
 *
 * Category: Financial
 *
 * @param {*} effect_rate The effective interest rate.
 * @param {*} npery The number of compounding periods per year.
 * @returns
 */
function NOMINAL(effect_rate, npery) {
  effect_rate = parseNumber(effect_rate);
  npery = parseNumber(npery);

  if (anyIsError(effect_rate, npery)) {
    return value
  }

  // Return error if rate <=0 or periods < 1
  if (effect_rate <= 0 || npery < 1) {
    return num
  }

  // Truncate periods if it is not an integer
  npery = parseInt(npery, 10);

  // Return nominal annual interest rate
  return (Math.pow(effect_rate + 1, 1 / npery) - 1) * npery
}

/**
 * Returns the number of periods for an investment.
 *
 * Category: Financial
 *
 * @param {*} rate The interest rate per period.
 * @param {*} pmt The payment made each period; it cannot change over the life of the annuity. Typically, pmt contains principal and interest but no other fees or taxes.
 * @param {*} pv The present value, or the lump-sum amount that a series of future payments is worth right now.
 * @param {*} fv Optional. The future value, or a cash balance you want to attain after the last payment is made. If fv is omitted, it is assumed to be 0 (the future value of a loan, for example, is 0).
 * @param {*} type Optional. The number 0 or 1 and indicates when payments are due.
 * @returns
 */
function NPER(rate, pmt, pv, fv, type) {
  type = type === undefined ? 0 : type;
  fv = fv === undefined ? 0 : fv;

  rate = parseNumber(rate);
  pmt = parseNumber(pmt);
  pv = parseNumber(pv);
  fv = parseNumber(fv);
  type = parseNumber(type);

  if (anyIsError(rate, pmt, pv, fv, type)) {
    return value
  }

  if (rate === 0) {
    return -(pv + fv) / pmt
  } else {
    const num = pmt * (1 + rate * type) - fv * rate;
    const den = pv * rate + pmt * (1 + rate * type);

    return Math.log(num / den) / Math.log(1 + rate)
  }
}

/**
 * Returns the net present value of an investment based on a series of periodic cash flows and a discount rate.
 *
 * Category: Financial
 *
 * @param {*} rate The rate of discount over the length of one period.
 * @param {*} args value1, value2, ... Value1 is required, subsequent values are optional. 1 to 254 arguments representing the payments and income.
 - value1, value2, ... must be equally spaced in time and occur at the end of each period.
 - NPV uses the order of value1, value2, ... to interpret the order of cash flows. Be sure to enter your payment and income values in the correct sequence.
 - Arguments that are empty values, logical values, or text representations of numbers, error values, or text that cannot be translated into numbers are ignored.
 - If an argument is an array or reference, only numbers in that array or reference are counted. Empty values, logical values, text, or error values in the array or reference are ignored.
 * @returns
 */
function NPV() {
  const args = parseNumberArray(flatten(arguments));

  if (args instanceof Error) {
    return args
  }

  // Lookup rate
  const rate = args[0];

  // Initialize net present value
  let value = 0;

  // Loop on all values
  for (let j = 1; j < args.length; j++) {
    value += args[j] / Math.pow(1 + rate, j);
  }

  // Return net present value
  return value
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the price per $100 face value of a security with an odd first period.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} issue The security's issue date.
 * @param {*} first_coupon The security's first coupon date.
 * @param {*} rate The security's interest rate.
 * @param {*} yld The security's annual yield.
 * @param {*} redemption The security's redemption value per $100 face value.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function ODDFPRICE() {
  throw new Error('ODDFPRICE is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the yield of a security with an odd first period.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} issue The security's issue date.
 * @param {*} first_coupon The security's first coupon date.
 * @param {*} rate The security's interest rate.
 * @param {*} pr The security's price.
 * @param {*} redemption The security's redemption value per $100 face value.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function ODDFYIELD() {
  throw new Error('ODDFYIELD is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the price per $100 face value of a security with an odd last period.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} last_interest The security's last coupon date.
 * @param {*} rate The security's interest rate.
 * @param {*} yld The security's annual yield.
 * @param {*} redemption The security's redemption value per $100 face value.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function ODDLPRICE() {
  throw new Error('ODDLPRICE is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the yield of a security with an odd last period.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} last_interest The security's last coupon date.
 * @param {*} rate The security's interest rate
 * @param {*} pr The security's price.
 * @param {*} redemption The security's redemption value per $100 face value.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function ODDLYIELD() {
  throw new Error('ODDLYIELD is not implemented')
}

/**
 * Returns the number of periods required by an investment to reach a specified value.
 *
 * Category: Financial
 *
 * @param {*} rate Rate is the interest rate per period.
 * @param {*} pv Pv is the present value of the investment.
 * @param {*} fv Fv is the desired future value of the investment.
 * @returns
 */
function PDURATION(rate, pv, fv) {
  rate = parseNumber(rate);
  pv = parseNumber(pv);
  fv = parseNumber(fv);

  if (anyIsError(rate, pv, fv)) {
    return value
  }

  // Return error if rate <=0
  if (rate <= 0) {
    return num
  }

  // Return number of periods
  return (Math.log(fv) - Math.log(pv)) / Math.log(1 + rate)
}

/**
 * Returns the periodic payment for an annuity.
 *
 * Category: Financial
 *
 * @param {*} rate The interest rate for the loan.
 * @param {*} nper The total number of payments for the loan.
 * @param {*} pv The present value, or the total amount that a series of future payments is worth now; also known as the principal.
 * @param {*} fv Optional. The future value, or a cash balance you want to attain after the last payment is made. If fv is omitted, it is assumed to be 0 (zero), that is, the future value of a loan is 0.
 * @param {*} type Optional. The number 0 (zero) or 1 and indicates when payments are due.
 * @returns
 */
function PMT(rate, nper, pv, fv, type) {
  // Credits: algorithm inspired by Apache OpenOffice
  fv = fv || 0;
  type = type || 0;

  rate = parseNumber(rate);
  nper = parseNumber(nper);
  pv = parseNumber(pv);
  fv = parseNumber(fv);
  type = parseNumber(type);

  if (anyIsError(rate, nper, pv, fv, type)) {
    return value
  }

  // Return payment
  let result;

  if (rate === 0) {
    result = (pv + fv) / nper;
  } else {
    const term = Math.pow(1 + rate, nper);

    result =
      type === 1
        ? ((fv * rate) / (term - 1) + (pv * rate) / (1 - 1 / term)) / (1 + rate)
        : (fv * rate) / (term - 1) + (pv * rate) / (1 - 1 / term);
  }

  return -result
}

/**
 * Returns the payment on the principal for an investment for a given period.
 *
 * Category: Financial
 *
 * @param {*} rate The interest rate per period.
 * @param {*} per Specifies the period and must be in the range 1 to nper.
 * @param {*} nper The total number of payment periods in an annuity.
 * @param {*} pv The present value — the total amount that a series of future payments is worth now.
 * @param {*} fv Optional. The future value, or a cash balance you want to attain after the last payment is made. If fv is omitted, it is assumed to be 0 (zero), that is, the future value of a loan is 0.
 * @param {*} type Optional. The number 0 or 1 and indicates when payments are due.
 * @returns
 */
function PPMT(rate, per, nper, pv, fv, type) {
  fv = fv || 0;
  type = type || 0;

  rate = parseNumber(rate);
  nper = parseNumber(nper);
  pv = parseNumber(pv);
  fv = parseNumber(fv);
  type = parseNumber(type);

  if (anyIsError(rate, nper, pv, fv, type)) {
    return value
  }

  return PMT(rate, nper, pv, fv, type) - IPMT(rate, per, nper, pv, fv, type)
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the price per $100 face value of a security that pays periodic interest.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} rate The security's annual coupon rate.
 * @param {*} yld The security's annual yield.
 * @param {*} redemption The security's redemption value per $100 face value.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function PRICE() {
  throw new Error('PRICE is not implemented')
}

/**
 * Returns the price per $100 face value of a discounted security.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} discount The security's discount rate.
 * @param {*} redemption The security's redemption value per $100 face value.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function PRICEDISC(settlement, maturity, discount, redemption, basis) {
  settlement = parseDate(settlement);
  maturity = parseDate(maturity);
  discount = parseNumber(discount);
  redemption = parseNumber(redemption);
  basis = parseNumber(basis);

  basis = basis || 0;

  if (anyIsError(settlement, maturity, discount, redemption, basis)) {
    return value
  }

  if (discount <= 0 || redemption <= 0) {
    return num
  }

  if (settlement >= maturity) {
    return value
  }

  let basisVal, diff;
  switch (basis) {
    case 0:
      basisVal = 360;
      diff = DAYS360(settlement, maturity, false);
      break
    case 1:
      basisVal = 365;
      diff = DATEDIF(settlement, maturity, 'D');
      break
    case 2:
      basisVal = 360;
      diff = DATEDIF(settlement, maturity, 'D');
      break
    case 3:
      basisVal = 365;
      diff = DATEDIF(settlement, maturity, 'D');
      break
    case 4:
      basisVal = 360;
      diff = DAYS360(settlement, maturity, true);
      break
    default:
      return num
  }

  return redemption - (discount * redemption * diff) / basisVal
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the price per $100 face value of a security that pays interest at maturity.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} issue The security's issue date, expressed as a serial date number.
 * @param {*} rate The security's interest rate at date of issue.
 * @param {*} yld The security's annual yield.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function PRICEMAT() {
  throw new Error('PRICEMAT is not implemented')
}

/**
 * Returns the present value of an investment.
 *
 * Category: Financial
 *
 * @param {*} rate The interest rate per period. For example, if you obtain an automobile loan at a 10 percent annual interest rate and make monthly payments, your interest rate per month is 10%/12, or 0.83%. You would enter 10%/12, or 0.83%, or 0.0083, into the formula as the rate.
 * @param {*} nper The total number of payment periods in an annuity. For example, if you get a four-year car loan and make monthly payments, your loan has 4*12 (or 48) periods. You would enter 48 into the formula for nper.
 * @param {*} pmt The payment made each period and cannot change over the life of the annuity. Typically, pmt includes principal and interest but no other fees or taxes. For example, the monthly payments on a $10,000, four-year car loan at 12 percent are $263.33. You would enter -263.33 into the formula as the pmt. If pmt is omitted, you must include the fv argument.
 * @param {*} fv Optional. The future value, or a cash balance you want to attain after the last payment is made. If fv is omitted, it is assumed to be 0 (the future value of a loan, for example, is 0). For example, if you want to save $50,000 to pay for a special project in 18 years, then $50,000 is the future value. You could then make a conservative guess at an interest rate and determine how much you must save each month. If fv is omitted, you must include the pmt argument.
 * @param {*} type Optional. The number 0 or 1 and indicates when payments are due.
 * @returns
 */
function PV(rate, per, pmt, fv, type) {
  fv = fv || 0;
  type = type || 0;

  rate = parseNumber(rate);
  per = parseNumber(per);
  pmt = parseNumber(pmt);
  fv = parseNumber(fv);
  type = parseNumber(type);

  if (anyIsError(rate, per, pmt, fv, type)) {
    return value
  }

  // Return present value
  return rate === 0
    ? -pmt * per - fv
    : (((1 - Math.pow(1 + rate, per)) / rate) * pmt * (1 + rate * type) - fv) / Math.pow(1 + rate, per)
}

/**
 * Returns the interest rate per period of an annuity.
 *
 * Category: Financial
 *
 * @param {*} nper The total number of payment periods in an annuity.
 * @param {*} pmt The payment made each period and cannot change over the life of the annuity. Typically, pmt includes principal and interest but no other fees or taxes. If pmt is omitted, you must include the fv argument.
 * @param {*} pv The present value — the total amount that a series of future payments is worth now.
 * @param {*} fv Optional. The future value, or a cash balance you want to attain after the last payment is made. If fv is omitted, it is assumed to be 0 (the future value of a loan, for example, is 0). If fv is omitted, you must include the pmt argument.
 * @param {*} type Optional. The number 0 or 1 and indicates when payments are due.
 * @param {*} guess Optional. Your guess for what the rate will be. If you omit guess, it is assumed to be 10 percent. If RATE does not converge, try different values for guess. RATE usually converges if guess is between 0 and 1.
 - If you omit guess, it is assumed to be 10 percent.
 - If RATE does not converge, try different values for guess. RATE usually converges if guess is between 0 and 1.
 * @returns
 */
function RATE(nper, pmt, pv, fv, type, guess) {
  guess = guess === undefined ? 0.01 : guess;
  fv = fv === undefined ? 0 : fv;
  type = type === undefined ? 0 : type;

  nper = parseNumber(nper);
  pmt = parseNumber(pmt);
  pv = parseNumber(pv);
  fv = parseNumber(fv);
  type = parseNumber(type);
  guess = parseNumber(guess);

  if (anyIsError(nper, pmt, pv, fv, type, guess)) {
    return value
  }

  const epsMax = 1e-10;
  const iterMax = 20;
  let rate = guess;

  type = type ? 1 : 0;

  for (let i = 0; i < iterMax; i++) {
    if (rate <= -1) {
      return num
    }

    let y, f;

    if (Math.abs(rate) < epsMax) {
      y = pv * (1 + nper * rate) + pmt * (1 + rate * type) * nper + fv;
    } else {
      f = Math.pow(1 + rate, nper);
      y = pv * f + pmt * (1 / rate + type) * (f - 1) + fv;
    }

    if (Math.abs(y) < epsMax) {
      return rate
    }

    let dy;

    if (Math.abs(rate) < epsMax) {
      dy = pv * nper + pmt * type * nper;
    } else {
      f = Math.pow(1 + rate, nper);
      const df = nper * Math.pow(1 + rate, nper - 1);
      dy = pv * df + pmt * (1 / rate + type) * df + pmt * (-1 / (rate * rate)) * (f - 1);
    }

    rate -= y / dy;
  }

  return rate
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the amount received at maturity for a fully invested security.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} investment The amount invested in the security.
 * @param {*} discount The security's discount rate.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function RECEIVED() {
  throw new Error('RECEIVED is not implemented')
}

/**
 * Returns an equivalent interest rate for the growth of an investment.
 *
 * Category: Financial
 *
 * @param {*} nper Nper is the number of periods for the investment.
 * @param {*} pv Pv is the present value of the investment.
 * @param {*} fv Fv is the future value of the investment.
 * @returns
 */
function RRI(nper, pv, fv) {
  nper = parseNumber(nper);
  pv = parseNumber(pv);
  fv = parseNumber(fv);

  if (anyIsError(nper, pv, fv)) {
    return value
  }

  // Return error if nper or present is equal to 0 (zero)
  if (nper === 0 || pv === 0) {
    return num
  }

  // Return equivalent interest rate
  return Math.pow(fv / pv, 1 / nper) - 1
}

/**
 * Returns the straight-line depreciation of an asset for one period.
 *
 * Category: Financial
 *
 * @param {*} cost The initial cost of the asset.
 * @param {*} salvage The value at the end of the depreciation (sometimes called the salvage value of the asset).
 * @param {*} life The number of periods over which the asset is depreciated (sometimes called the useful life of the asset).
 * @returns
 */
function SLN(cost, salvage, life) {
  cost = parseNumber(cost);
  salvage = parseNumber(salvage);
  life = parseNumber(life);

  if (anyIsError(cost, salvage, life)) {
    return value
  }

  // Return error if life equal to 0 (zero)
  if (life === 0) {
    return num
  }

  // Return straight-line depreciation
  return (cost - salvage) / life
}

/**
 * Returns the sum-of-years' digits depreciation of an asset for a specified period.
 *
 * Category: Financial
 *
 * @param {*} cost The initial cost of the asset.
 * @param {*} salvage The value at the end of the depreciation (sometimes called the salvage value of the asset).
 * @param {*} life The number of periods over which the asset is depreciated (sometimes called the useful life of the asset).
 * @param {*} per The period and must use the same units as life.
 * @returns
 */
function SYD(cost, salvage, life, per) {
  // Return error if any of the parameters is not a number
  cost = parseNumber(cost);
  salvage = parseNumber(salvage);
  life = parseNumber(life);
  per = parseNumber(per);

  if (anyIsError(cost, salvage, life, per)) {
    return value
  }

  // Return error if life equal to 0 (zero)
  if (life === 0) {
    return num
  }

  // Return error if period is lower than 1 or greater than life
  if (per < 1 || per > life) {
    return num
  }

  // Truncate period if it is not an integer
  per = parseInt(per, 10);

  // Return straight-line depreciation
  return ((cost - salvage) * (life - per + 1) * 2) / (life * (life + 1))
}

/**
 * Returns the bond-equivalent yield for a Treasury bill.
 *
 * Category: Financial
 *
 * @param {*} settlement The Treasury bill's settlement date. The security settlement date is the date after the issue date when the Treasury bill is traded to the buyer.
 * @param {*} maturity The Treasury bill's maturity date. The maturity date is the date when the Treasury bill expires.
 * @param {*} discount The Treasury bill's discount rate.
 * @returns
 */
function TBILLEQ(settlement, maturity, discount) {
  settlement = parseDate(settlement);
  maturity = parseDate(maturity);
  discount = parseNumber(discount);

  if (anyIsError(settlement, maturity, discount)) {
    return value
  }

  // Return error if discount is lower than or equal to zero
  if (discount <= 0) {
    return num
  }

  // Return error if settlement is greater than maturity
  if (settlement > maturity) {
    return num
  }

  // Return error if maturity is more than one year after settlement
  if (maturity - settlement > 365 * 24 * 60 * 60 * 1000) {
    return num
  }

  // Return bond-equivalent yield
  return (365 * discount) / (360 - discount * DAYS360(settlement, maturity, false))
}

/**
 * Returns the price per $100 face value for a Treasury bill.
 *
 * Category: Financial
 *
 * @param {*} settlement The Treasury bill's settlement date. The security settlement date is the date after the issue date when the Treasury bill is traded to the buyer.
 * @param {*} maturity The Treasury bill's maturity date. The maturity date is the date when the Treasury bill expires.
 * @param {*} discount The Treasury bill's discount rate.
 * @returns
 */
function TBILLPRICE(settlement, maturity, discount) {
  settlement = parseDate(settlement);
  maturity = parseDate(maturity);
  discount = parseNumber(discount);

  if (anyIsError(settlement, maturity, discount)) {
    return value
  }

  // Return error if discount is lower than or equal to zero
  if (discount <= 0) {
    return num
  }

  // Return error if settlement is greater than maturity
  if (settlement > maturity) {
    return num
  }

  // Return error if maturity is more than one year after settlement
  if (maturity - settlement > 365 * 24 * 60 * 60 * 1000) {
    return num
  }

  // Return bond-equivalent yield
  return 100 * (1 - (discount * DAYS360(settlement, maturity, false)) / 360)
}

/**
 * Returns the yield for a Treasury bill.
 *
 * Category: Financial
 *
 * @param {*} settlement The Treasury bill's settlement date. The security settlement date is the date after the issue date when the Treasury bill is traded to the buyer.
 * @param {*} maturity The Treasury bill's maturity date. The maturity date is the date when the Treasury bill expires.
 * @param {*} pr The Treasury bill's price per $100 face value.
 * @returns
 */
function TBILLYIELD(settlement, maturity, pr) {
  settlement = parseDate(settlement);
  maturity = parseDate(maturity);
  pr = parseNumber(pr);

  if (anyIsError(settlement, maturity, pr)) {
    return value
  }

  // Return error if price is lower than or equal to zero
  if (pr <= 0) {
    return num
  }

  // Return error if settlement is greater than maturity
  if (settlement > maturity) {
    return num
  }

  // Return error if maturity is more than one year after settlement
  if (maturity - settlement > 365 * 24 * 60 * 60 * 1000) {
    return num
  }

  // Return bond-equivalent yield
  return ((100 - pr) * 360) / (pr * DAYS360(settlement, maturity, false))
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the depreciation of an asset for a specified or partial period by using a declining balance method.
 *
 * Category: Financial
 *
 * @param {*} cost The initial cost of the asset.
 * @param {*} salvage The value at the end of the depreciation (sometimes called the salvage value of the asset). This value can be 0.
 * @param {*} life The number of periods over which the asset is depreciated (sometimes called the useful life of the asset).
 * @param {*} start_period The starting period for which you want to calculate the depreciation. Start_period must use the same units as life.
 * @param {*} end_period The ending period for which you want to calculate the depreciation. End_period must use the same units as life.
 * @param {*} factor Optional. The rate at which the balance declines. If factor is omitted, it is assumed to be 2 (the double-declining balance method). Change factor if you do not want to use the double-declining balance method. For a description of the double-declining balance method, see DDB.
 * @param {*} no_switch Optional. A logical value specifying whether to switch to straight-line depreciation when depreciation is greater than the declining balance calculation.
 - If no_switch is TRUE, Microsoft Excel does not switch to straight-line depreciation even when the depreciation is greater than the declining balance calculation.
 - If no_switch is FALSE or omitted, Excel switches to straight-line depreciation when depreciation is greater than the declining balance calculation.
 * @returns
 */
function VDB() {
  throw new Error('VDB is not implemented')
}

/**
 * Returns the internal rate of return for a schedule of cash flows that is not necessarily periodic.
 *
 * Category: Financial
 *
 * @param {*} values A series of cash flows that corresponds to a schedule of payments in dates. The first payment is optional and corresponds to a cost or payment that occurs at the beginning of the investment. If the first value is a cost or payment, it must be a negative value. All succeeding payments are discounted based on a 365-day year. The series of values must contain at least one positive and one negative value.
 * @param {*} dates A schedule of payment dates that corresponds to the cash flow payments. Dates may occur in any order. Dates should be entered by using the DATE function, or as results of other formulas or functions. For example, use DATE(2008,5,23) for the 23rd day of May, 2008. Problems can occur if dates are entered as text. .
 * @param {*} guess Optional. A number that you guess is close to the result of XIRR.
 * @returns
 */
function XIRR(values, dates, guess) {
  // Credits: algorithm inspired by Apache OpenOffice
  values = parseNumberArray(flatten(values));
  dates = parseDateArray(flatten(dates));
  guess = parseNumber(guess);

  if (anyIsError(values, dates, guess)) {
    return value
  }

  // Calculates the resulting amount
  const irrResult = (values, dates, rate) => {
    const r = rate + 1;
    let result = values[0];

    for (let i = 1; i < values.length; i++) {
      result += values[i] / Math.pow(r, DAYS(dates[i], dates[0]) / 365);
    }

    return result
  };

  // Calculates the first derivation
  const irrResultDeriv = (values, dates, rate) => {
    const r = rate + 1;
    let result = 0;

    for (let i = 1; i < values.length; i++) {
      const frac = DAYS(dates[i], dates[0]) / 365;
      result -= (frac * values[i]) / Math.pow(r, frac + 1);
    }

    return result
  };

  // Check that values contains at least one positive value and one negative value
  let positive = false;
  let negative = false;

  for (let i = 0; i < values.length; i++) {
    if (values[i] > 0) {
      positive = true;
    }

    if (values[i] < 0) {
      negative = true;
    }
  }

  // Return error if values does not contain at least one positive value and one negative value
  if (!positive || !negative) {
    return num
  }

  // Initialize guess and resultRate
  guess = guess || 0.1;
  let resultRate = guess;

  // Set maximum epsilon for end of iteration
  const epsMax = 1e-10;

  // Implement Newton's method
  let newRate, epsRate, resultValue;
  let contLoop = true;

  do {
    resultValue = irrResult(values, dates, resultRate);
    newRate = resultRate - resultValue / irrResultDeriv(values, dates, resultRate);
    epsRate = Math.abs(newRate - resultRate);
    resultRate = newRate;
    contLoop = epsRate > epsMax && Math.abs(resultValue) > epsMax;
  } while (contLoop)

  // Return internal rate of return
  return resultRate
}

/**
 * Returns the net present value for a schedule of cash flows that is not necessarily periodic.
 *
 * Category: Financial
 *
 * @param {*} rate The discount rate to apply to the cash flows.
 * @param {*} values A series of cash flows that corresponds to a schedule of payments in dates. The first payment is optional and corresponds to a cost or payment that occurs at the beginning of the investment. If the first value is a cost or payment, it must be a negative value. All succeeding payments are discounted based on a 365-day year. The series of values must contain at least one positive value and one negative value.
 * @param {*} dates A schedule of payment dates that corresponds to the cash flow payments. The first payment date indicates the beginning of the schedule of payments. All other dates must be later than this date, but they may occur in any order.
 * @returns
 */
function XNPV(rate, values, dates) {
  rate = parseNumber(rate);
  values = parseNumberArray(flatten(values));
  dates = parseDateArray(flatten(dates));

  if (anyIsError(rate, values, dates)) {
    return value
  }

  let result = 0;

  for (let i = 0; i < values.length; i++) {
    result += values[i] / Math.pow(1 + rate, DAYS(dates[i], dates[0]) / 365);
  }

  return result
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the yield on a security that pays periodic interest.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} rate The security's annual coupon rate.
 * @param {*} pr The security's price per $100 face value.
 * @param {*} redemption The security's redemption value per $100 face value.
 * @param {*} frequency The number of coupon payments per year. For annual payments, frequency = 1; for semiannual, frequency = 2; for quarterly, frequency = 4.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function YIELD() {
  throw new Error('YIELD is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the annual yield for a discounted security; for example, a Treasury bill.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} pr The security's price per $100 face value.
 * @param {*} redemption The security's redemption value per $100 face value.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function YIELDDISC() {
  throw new Error('YIELDDISC is not implemented')
}

// TODO
/**
 * -- Not implemented --
 *
 * Returns the annual yield of a security that pays interest at maturity.
 *
 * Category: Financial
 *
 * @param {*} settlement The security's settlement date. The security settlement date is the date after the issue date when the security is traded to the buyer.
 * @param {*} maturity The security's maturity date. The maturity date is the date when the security expires.
 * @param {*} issue The security's issue date, expressed as a serial date number.
 * @param {*} rate The security's interest rate at date of issue.
 * @param {*} pr The security's price per $100 face value.
 * @param {*} basis Optional. The type of day count basis to use.
 * @returns
 */
function YIELDMAT() {
  throw new Error('YIELDMAT is not implemented')
}

/**
 * Returns TRUE if all of its arguments are TRUE.
 *
 * Category: Logical
 *
 * @returns
 */
function AND() {
  const args = flatten(arguments);
  let result = value;

  for (let i = 0; i < args.length; i++) {
    if (args[i] instanceof Error) {
      return args[i]
    }

    if (args[i] === undefined || args[i] === null || typeof args[i] === 'string') {
      continue
    }

    if (result === value) {
      result = true;
    }

    if (!args[i]) {
      result = false;
    }
  }

  return result
}

/**
 * Returns the logical value FALSE.
 *
 * Category: Logical
 *
 * @returns
 */
function FALSE() {
  return false
}

/**
 * Specifies a logical test to perform.
 *
 * Category: Logical
 *
 * @param {*} logical_test
 * @param {*} value_if_true
 * @param {*} value_if_false
 *
 * @returns
 */
function IF(logical_test, value_if_true, value_if_false) {
  if (logical_test instanceof Error) {
    return logical_test
  }

  value_if_true = arguments.length >= 2 ? value_if_true : true;

  if (value_if_true === undefined || value_if_true === null) {
    value_if_true = 0;
  }

  value_if_false = arguments.length === 3 ? value_if_false : false;

  if (value_if_false === undefined || value_if_false === null) {
    value_if_false = 0;
  }

  return logical_test ? value_if_true : value_if_false
}

/**
 * Checks whether one or more conditions are met and returns a value that corresponds to the first TRUE condition.
 *
 * Category: Logical
 *
 * @returns
 */
function IFS() {
  for (let i = 0; i < arguments.length / 2; i++) {
    if (arguments[i * 2]) {
      return arguments[i * 2 + 1]
    }
  }

  return na
}

/**
 * Returns a value you specify if a formula evaluates to an error; otherwise, returns the result of the formula.
 *
 * Category: Logical
 *
 * @param {*} value The argument that is checked for an error.
 * @param {*} value_if_error The value to return if the formula evaluates to an error. The following error types are evaluated: #N/A, #VALUE!, #REF!, #DIV/0!, #NUM!, #NAME?, or #NULL!.
 * @returns
 */
function IFERROR(value, value_if_error) {
  if (ISERROR(value)) {
    return value_if_error
  }

  return value
}

/**
 * Returns the value you specify if the expression resolves to #N/A, otherwise returns the result of the expression.
 *
 * Category: Logical
 *
 * @returns
 */
function IFNA(value, value_if_na) {
  return value === na ? value_if_na : value
}

/**
 * Reverses the logic of its argument.
 *
 * Category: Logical
 *
 * @returns
 */
function NOT(logical) {
  if (typeof logical === 'string') {
    return value
  }

  if (logical instanceof Error) {
    return logical
  }

  return !logical
}

/**
 * Returns TRUE if any argument is TRUE.
 *
 * Category: Logical
 *
 * @returns
 */
function OR() {
  const args = flatten(arguments);
  let result = value;

  for (let i = 0; i < args.length; i++) {
    if (args[i] instanceof Error) {
      return args[i]
    }

    if (args[i] === undefined || args[i] === null || typeof args[i] === 'string') {
      continue
    }

    if (result === value) {
      result = false;
    }

    if (args[i]) {
      result = true;
    }
  }

  return result
}

/**
 * Returns the logical value TRUE.
 *
 * Category: Logical
 *
 * @returns
 */
function TRUE() {
  return true
}

/**
 * Returns a logical exclusive OR of all arguments.
 *
 * Category: Logical
 *
 * @param {*} args logical1, logical2,… Logical 1 is required, subsequent logical values are optional. 1 to 254 conditions you want to test that can be either TRUE or FALSE, and can be logical values, arrays, or references.
 * @returns
 */
function XOR() {
  const args = flatten(arguments);
  let result = value;

  for (let i = 0; i < args.length; i++) {
    if (args[i] instanceof Error) {
      return args[i]
    }

    if (args[i] === undefined || args[i] === null || typeof args[i] === 'string') {
      continue
    }

    if (result === value) {
      result = 0;
    }

    if (args[i]) {
      result++;
    }
  }

  if (result === value) {
    return result
  }

  return !!(Math.floor(Math.abs(result)) & 1)
}

/**
 * Evaluates an expression against a list of values and returns the result corresponding to the first matching value. If there is no match, an optional default value may be returned.
 *
 * Category: Logical
 *
 * @returns
 */
function SWITCH() {
  let result;

  if (arguments.length > 0) {
    const targetValue = arguments[0];
    const argc = arguments.length - 1;
    const switchCount = Math.floor(argc / 2);
    let switchSatisfied = false;
    const hasDefaultClause = argc % 2 !== 0;
    const defaultClause = argc % 2 === 0 ? null : arguments[arguments.length - 1];

    if (switchCount) {
      for (let index = 0; index < switchCount; index++) {
        if (targetValue === arguments[index * 2 + 1]) {
          result = arguments[index * 2 + 2];
          switchSatisfied = true;
          break
        }
      }
    }

    if (!switchSatisfied) {
      result = hasDefaultClause ? defaultClause : na;
    }
  } else {
    result = value;
  }

  return result
}

const utils = { errors, symbols };



;// CONCATENATED MODULE: ./src/formula.js
/**
 * Jspreadsheet Extensions
 * Extension: Formula Basic
 * License: This is a free software MIT
 *
 * https://jspreadsheet.com
 */



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

    for (let i = 0; i < Object.keys(esm_namespaceObject).length; i++) {
        let method = Object.keys(esm_namespaceObject)[i]
        let keys = []
        let values
        if (typeof esm_namespaceObject[method] == 'object') {
            keys = Object.keys(esm_namespaceObject[method])
            values = Object.values(esm_namespaceObject[method])
            for (let a = 0; a < values.length; a++) {
                if (typeof values[a] == 'object') {
                    let subMethod = keys[a]
                    if (esm_namespaceObject[method][subMethod]) {
                        keys = [
                            ...keys,
                            ...Object.keys(esm_namespaceObject[method][subMethod]).map((a) => subMethod + '.' + a)
                        ] // Line too heavy, need refactor
                        keys.splice(keys.indexOf(subMethod), 1)
                    }
                }
            }
        }

        if (keys.length < 1) {
            scope[method] = esm_namespaceObject[method]
        } else {
            for (let j = 0; j < keys.length; j++) {
                if (typeof getValue(esm_namespaceObject[method], keys[j]) == 'function') {
                    scope[method] = getValue(esm_namespaceObject[method], keys[j])
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

    /**
     * Get letter based on a number
     * @param {number} i
     * @return {string}
     */
    var getColumnName = function (i) {
        var letter = ''
        if (i > 701) {
            letter += String.fromCharCode(64 + parseInt(i / 676))
            letter += String.fromCharCode(64 + parseInt((i % 676) / 26))
        } else if (i > 25) {
            letter += String.fromCharCode(64 + parseInt(i / 26))
        }
        letter += String.fromCharCode(65 + (i % 26))

        return letter
    }

    /**
     * Get column name from coords
     */
    F.getColumnNameFromCoords = function (x, y) {
        return getColumnName(parseInt(x)) + (parseInt(y) + 1)
    }

    F.getCoordsFromColumnName = function (columnName) {
        // Get the letters
        var t = /^[a-zA-Z]+/.exec(columnName)

        if (t) {
            // Base 26 calculation
            var code = 0
            for (var i = 0; i < t[0].length; i++) {
                code += parseInt(t[0].charCodeAt(i) - 64) * Math.pow(26, t[0].length - 1 - i)
            }
            code--
            // Make sure jspreadsheet starts on zero
            if (code < 0) {
                code = 0
            }

            // Number
            var number = parseInt(/[0-9]+$/.exec(columnName)) || null
            if (number > 0) {
                number--
            }

            return [code, number]
        }
    }

    F.getRangeFromTokens = function (tokens) {
        tokens = tokens.filter(function (v) {
            return v != '#REF!'
        })

        var d = ''
        var t = ''
        for (var i = 0; i < tokens.length; i++) {
            if (tokens[i].indexOf('.') >= 0) {
                d = '.'
            } else if (tokens[i].indexOf('!') >= 0) {
                d = '!'
            }
            if (d) {
                t = tokens[i].split(d)
                tokens[i] = t[1]
                t = t[0] + d
            }
        }

        tokens.sort(function (a, b) {
            var t1 = Helpers.getCoordsFromColumnName(a)
            var t2 = Helpers.getCoordsFromColumnName(b)
            if (t1[1] > t2[1]) {
                return 1
            } else if (t1[1] < t2[1]) {
                return -1
            } else {
                if (t1[0] > t2[0]) {
                    return 1
                } else if (t1[0] < t2[0]) {
                    return -1
                } else {
                    return 0
                }
            }
        })

        if (!tokens.length) {
            return '#REF!'
        } else {
            return t + (tokens[0] + ':' + tokens[tokens.length - 1])
        }
    }

    F.getTokensFromRange = function (range) {
        if (range.indexOf('.') > 0) {
            var t = range.split('.')
            range = t[1]
            t = t[0] + '.'
        } else if (range.indexOf('!') > 0) {
            var t = range.split('!')
            range = t[1]
            t = t[0] + '!'
        } else {
            var t = ''
        }

        var range = range.split(':')
        var e1 = F.getCoordsFromColumnName(range[0])
        var e2 = F.getCoordsFromColumnName(range[1])

        if (e1[0] <= e2[0]) {
            var x1 = e1[0]
            var x2 = e2[0]
        } else {
            var x1 = e2[0]
            var x2 = e1[0]
        }

        if (e1[1] === null && e2[1] == null) {
            var y1 = null
            var y2 = null

            var k = Object.keys(vars)
            for (var i = 0; i < k.length; i++) {
                var tmp = F.getCoordsFromColumnName(k[i])
                if (tmp[0] === e1[0]) {
                    if (y1 === null || tmp[1] < y1) {
                        y1 = tmp[1]
                    }
                }
                if (tmp[0] === e2[0]) {
                    if (y2 === null || tmp[1] > y2) {
                        y2 = tmp[1]
                    }
                }
            }
        } else {
            if (e1[1] <= e2[1]) {
                var y1 = e1[1]
                var y2 = e2[1]
            } else {
                var y1 = e2[1]
                var y2 = e1[1]
            }
        }

        var f = []
        for (var j = y1; j <= y2; j++) {
            var line = []
            for (var i = x1; i <= x2; i++) {
                line.push(t + F.getColumnNameFromCoords(i, j))
            }
            f.push(line)
        }

        return f
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

let formula = Formula(undefined)
/* harmony default export */ var src_formula = (formula);
}();
__webpack_exports__ = __webpack_exports__["default"];
/******/ 	return __webpack_exports__;
/******/ })()
;
});