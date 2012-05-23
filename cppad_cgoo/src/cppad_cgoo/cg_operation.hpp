#ifndef CPPAD_CG_OPERATION_INCLUDED
#define	CPPAD_CG_OPERATION_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    /**
     * Possible operations
     * 
     * \author Joao Leal
     */
    enum CGOpCode {
        CGAbsOp, //  abs(variable)
        CGAcosOp, // asin(variable)
        CGAddOp, //  a + b
        CGAsinOp, // asin(variable)
        CGAtanOp, // atan(variable)
        CGComOpLt, // result = left < right? trueCase: falseCase
        CGComOpLe, // result = left <= right? trueCase: falseCase
        CGComOpEq, // result = left == right? trueCase: falseCase
        CGComOpGe, // result = left >= right? trueCase: falseCase
        CGComOpGt, // result = left > right? trueCase: falseCase
        CGComOpNe, // result = left != right? trueCase: falseCase
        CGCoshOp, // cosh(variable)
        CGCosOp, //  cos(variable)
        CGDivOp, // a / b
        CGExpOp, //  exp(variable)
        CGInvOp, //                             independent variable
        CGLogOp, //  log(variable)
        CGMulOp, // a * b
        CGPowOp, //  pow(a,   b)
//        PriOp, //  PrintFor(text, parameter or variable, parameter or variable)
        CGSignOp, // result = (x > 0)? 1.0:((x == 0)? 0.0:-1)
        CGSinhOp, // sinh(variable)
        CGSinOp, //  sin(variable)
        CGSqrtOp, // sqrt(variable)
        CGSubOp, //  a - b
        CGTanhOp, //  tan(variable)
        CGTanOp, //  tan(variable)
        CGUnMinusOp // -(a)
    };

}

#endif
