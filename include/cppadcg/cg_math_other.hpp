#ifndef CPPAD_CG_MATH_OTHER_INCLUDED
#define	CPPAD_CG_MATH_OTHER_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    template <class Base>
    inline CG<Base> pow(const CG<Base> &x, const CG<Base> &y) {
        if (x.isParameter() && y.isParameter()) {
            return CG<Base > (pow(x.getParameterValue(), y.getParameterValue()));
        }

        CodeHandler<Base>* handler;
        if (y.isParameter()) {
            if (y.IdenticalZero()) {
                return CG<Base > (Base(1.0)); // does not consider that x could be infinity
            } else if (y.IdenticalOne()) {
                return CG<Base > (x);
            }
            handler = x.getCodeHandler();
        } else {
            handler = y.getCodeHandler();
        }

        return CG<Base>(*handler, new SourceCodeFragment<Base>(CGPowOp, x.argument(), y.argument()));
    }

    template <class Base>
    inline CG<Base> pow(const Base &x, const CG<Base> &y) {
        return pow(CG<Base > (x), y);
    }

    template <class Base>
    inline CG<Base> pow(const CG<Base> &x, const Base &y) {
        return pow(CG<Base > (x), y);
    }

    template <class Base>
    inline CG<Base> sign(const CG<Base> &x) {
        if (x.isParameter()) {
            if (x.getParameterValue() > Base(1.0)) {
                return CG<Base > (Base(1.0));
            } else if (x.getParameterValue() == Base(0.0)) {
                return CG<Base > (Base(0.0));
            } else {
                return CG<Base > (Base(-1.0));
            }
        }

        return CG<Base>(*x.getCodeHandler(), new SourceCodeFragment<Base>(CGSignOp, x.argument()));
    }

}

#endif

