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

#include <cppad_cgoo/cg_cg.hpp>

namespace CppAD {

    template<class Base>
    inline bool IdenticalZero(const CG<Base>& x) throw (CGException) {
        if (x.isParameter()) {
            return IdenticalZero(x.getParameterValue());
        } else {
            throw CGException("Invalid call to IdenticalZero: argument is not a parameter");
        }
    }

    template<class Base>
    inline bool IdenticalOne(const CG<Base>& x) throw (CGException) {
        if (x.isParameter()) {
            return IdenticalOne(x.getParameterValue());
        } else {
            throw CGException("Invalid call to IdenticalOne: argument is not a parameter");
        }
    }

    template<class Base>
    inline bool IdenticalEqualPar(const CG<Base>& x, const CG<Base>& y) throw (CGException) {
        if (x.isParameter()) {
            if (y.isParameter()) {
                return IdenticalEqualPar(x.getParameterValue(), y.getParameterValue());
            } else {
                throw CGException("Invalid call to IdenticalEqualPar: 2nd argument is not a parameter");
            }
        } else {
            throw CGException("Invalid call to IdenticalEqualPar: 1st argument is not a parameter");
        }
    }
}

#endif

