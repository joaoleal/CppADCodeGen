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

#define CPPAD_CG_CREATE_OPERATION(OpName)                                      \
    template<class Base>                                                       \
    inline CG<Base> OpName(const CG<Base>& var) {                              \
        if (var.isParameter()) {                                               \
            return CG<Base > (OpName(var.getParameterValue()));                \
        } else {                                                               \
            std::string operations = #OpName"(" + var.operations() + ")";      \
            return CG<Base>(*var.getCodeHandler(), operations, FUNCTION, &var);\
        }                                                                      \
    }

    CPPAD_CG_CREATE_OPERATION(abs);
    CPPAD_CG_CREATE_OPERATION(acos);
    CPPAD_CG_CREATE_OPERATION(asin);
    CPPAD_CG_CREATE_OPERATION(atan);
    CPPAD_CG_CREATE_OPERATION(cos);
    CPPAD_CG_CREATE_OPERATION(cosh);
    CPPAD_CG_CREATE_OPERATION(exp);
    CPPAD_CG_CREATE_OPERATION(log);
    CPPAD_CG_CREATE_OPERATION(sin);
    CPPAD_CG_CREATE_OPERATION(sinh);
    CPPAD_CG_CREATE_OPERATION(sqrt);
    CPPAD_CG_CREATE_OPERATION(tan);
    CPPAD_CG_CREATE_OPERATION(tanh);
}
#endif

