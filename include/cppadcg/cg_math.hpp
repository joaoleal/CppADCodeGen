#ifndef CPPAD_CG_MATH_INCLUDED
#define	CPPAD_CG_MATH_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

#define CPPAD_CG_CREATE_OPERATION(OpName, OpCode)                              \
    template<class Base>                                                       \
    inline CG<Base> OpName(const CG<Base>& var) {                              \
        if (var.isParameter()) {                                               \
            return CG<Base > (OpName(var.getParameterValue()));                \
        } else {                                                               \
            return CG<Base>(*var.getCodeHandler(), new SourceCodeFragment<Base>(OpCode, var.argument()));\
        }                                                                      \
    }

    CPPAD_CG_CREATE_OPERATION(abs, CGAbsOp)
    CPPAD_CG_CREATE_OPERATION(acos, CGAcosOp)
    CPPAD_CG_CREATE_OPERATION(asin, CGAsinOp)
    CPPAD_CG_CREATE_OPERATION(atan, CGAtanOp)
    CPPAD_CG_CREATE_OPERATION(cos, CGCosOp)
    CPPAD_CG_CREATE_OPERATION(cosh, CGCoshOp)
    CPPAD_CG_CREATE_OPERATION(exp, CGExpOp)
    CPPAD_CG_CREATE_OPERATION(log, CGLogOp)
    CPPAD_CG_CREATE_OPERATION(sin, CGSinOp)
    CPPAD_CG_CREATE_OPERATION(sinh, CGSinhOp)
    CPPAD_CG_CREATE_OPERATION(sqrt, CGSqrtOp)
    CPPAD_CG_CREATE_OPERATION(tan, CGTanOp)
    CPPAD_CG_CREATE_OPERATION(tanh, CGTanhOp)
}
#endif

