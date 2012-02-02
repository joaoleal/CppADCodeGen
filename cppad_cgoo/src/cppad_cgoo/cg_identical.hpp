#ifndef CPPAD_CG_IDENTICAL_INCLUDED
#define	CPPAD_CG_IDENTICAL_INCLUDED
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
    inline bool IdenticalPar(const CG<Base>& x) throw (CGException) {
        if (!x.isParameter()) {
            throw CGException("Invalid call to IdenticalPar(): argument is not a parameter");
        }
        return IdenticalPar(x.getParameterValue());
    }

    template<class Base>
    inline bool IdenticalZero(const CG<Base>& x) throw (CGException) {
        if (!x.isParameter()) {
            throw CGException("Invalid call to IdenticalZero(): argument is not a parameter");
        }
        return IdenticalZero(x.getParameterValue());
    }

    template<class Base>
    inline bool IdenticalOne(const CG<Base>& x) throw (CGException) {
        if (!x.isParameter()) {
            throw CGException("Invalid call to IdenticalOne(): argument is not a parameter");
        }
        return IdenticalOne(x.getParameterValue());
    }

    template<class Base>
    inline bool IdenticalEqualPar(const CG<Base>& x, const CG<Base>& y) {
        return x.isParameter() && y.isParameter() && IdenticalEqualPar(x.getParameterValue(), y.getParameterValue());
    }
}

#endif

