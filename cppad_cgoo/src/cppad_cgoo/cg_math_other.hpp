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
        CPPAD_CG_CHECK_CG(x);
        CPPAD_CG_CHECK_CG(y);

        if (x.isParameter() && y.isParameter()) {
            return CG<Base > (pow(x.getParameterValue(), y.getParameterValue()));
        }

        CodeHandler<Base>* handle;
        if (y.isParameter()) {
            if (y.IdenticalZero()) {
                return CG<Base > (Base(1.0)); // does not consider that x could be infinity
            } else if (y.IdenticalOne()) {
                return CG<Base > (x);
            }
            handle = x.getCodeHandler();
        } else {
            handle = y.getCodeHandler();
        }

        std::string operations = "pow(" + handle->operations(x) + ", " + handle->operations(y) + ")";
        CG<Base> result;
        result.makeTemporaryVariable(*handle, operations, FUNCTION);
        return result;
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
        CPPAD_CG_CHECK_CG(x);

        if (x.isParameter()) {
            if (x.getParameterValue() > Base(1.0)) {
                return CG<Base > (Base(1.0));
            } else if (x.getParameterValue() == Base(0.0)) {
                return CG<Base > (Base(0.0));
            } else {
                return CG<Base > (Base(-1.0));
            }
        }

        const CG<Base>* x1;
        CG<Base> y;
        if (x.isTemporaryVariable()) {
            // make it a variable
            y = x;
            x1 = &y;
        } else {
            x1 = &x;
        }

        std::string name = x1->createVariableName();
        CodeHandler<Base>* h = x1->getCodeHandler();
        std::string operations = "(" + name + " > " + h->baseToString(0.0) + "?"
                + h->baseToString(1.0) + ":("
                + name + " == " + h->baseToString(0.0) + "?"
                + h->baseToString(0.0) + ":"
                + h->baseToString(-1.0)
                + ")"
                ")";
        CG<Base> result;
        result.makeTemporaryVariable(*x1->getCodeHandler(), operations, FUNCTION);
        return result;
    }

}

#endif

