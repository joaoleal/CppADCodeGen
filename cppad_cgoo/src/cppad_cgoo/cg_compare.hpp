#ifndef CPPAD_CG_COMPARE_INCLUDED
#define	CPPAD_CG_COMPARE_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    template<class Base>
    inline bool operator ==(const CG<Base> &left, const CG<Base> &right) {
        if (left.isParameter() && right.isParameter()) {
            return left.getParameterValue() == right.getParameterValue();
        } else if (left.isParameter() || right.isParameter()) {
            return false;
        } else {
            return left.getVariableID() == right.getVariableID();
        }
    }

    template<class Base>
    inline bool operator !=(const CG<Base> &left, const CG<Base> &right) {
        if (left.isParameter() && right.isParameter()) {
            return left.getParameterValue() != right.getParameterValue();
        } else if (left.isParameter() || right.isParameter()) {
            return true;
        } else {
            return left.getVariableID() != right.getVariableID();
        }
    }

#define CPPAD_CG_OPERATOR(Op)                                                  \
    template<class Base>                                                       \
    inline bool operator Op(const CG<Base> &left, const CG<Base> &right) {     \
        if (left.isParameter() && right.isParameter()) {                       \
            return left.getParameterValue() Op right.getParameterValue();      \
        } else {                                                               \
            throw CGException("Cannot use the "#Op" comparison operator in non parameter variables");\
        }                                                                      \
    }

    CPPAD_CG_OPERATOR(>)
    CPPAD_CG_OPERATOR( >=)
    CPPAD_CG_OPERATOR(<)
    CPPAD_CG_OPERATOR( <=)
}

#endif

