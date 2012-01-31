/* 
 * File:   cg_ordered.hpp
 * Author: joao
 *
 * Created on 31 de Janeiro de 2012, 12:38
 */

#ifndef CPPAD_CG_ORDERED_INCLUDED
#define	CPPAD_CG_ORDERED_INCLUDED
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
    bool GreaterThanZero(const CG<Base> &x) {
        if (!x.isParameter()) {
            throw CGException("GreaterThanZero cannot be called for non-parameters");
        }

        return GreaterThanZero(x.getParameterValue());
    }

    template<class Base>
    bool GreaterThanOrZero(const CG<Base> &x) {
        if (!x.isParameter()) {
            throw CGException("GreaterThanOrZero cannot be called for non-parameters");
        }

        return GreaterThanOrZero(x.getParameterValue());
    }

    template<class Base>
    bool LessThanZero(const CG<Base> &x) {
        if (!x.isParameter()) {
            throw CGException("LessThanZero cannot be called for non-parameters");
        }

        return LessThanZero(x.getParameterValue());
    }

    template<class Base>
    bool LessThanOrZero(const CG<Base> &x) {
        if (!x.isParameter()) {
            throw CGException("LessThanOrZero cannot be called for non-parameters");
        }

        return LessThanOrZero(x.getParameterValue());
    }

    template<class Base>
    bool abs_geq(const CG<Base>& x, const CG<Base>& y) {
        if (!x.isParameter()) {
            throw CGException("abs_geq cannot be called for non-parameters (x)");
        } else if (!y.isParameter()) {
            throw CGException("abs_geq cannot be called for non-parameters (y)");
        }

        return abs_geq(x.getParameterValue(), y.getParameterValue());
    }

}

#endif

