#ifndef CPPAD_CG_MATH_OTHER_INCLUDED
#define CPPAD_CG_MATH_OTHER_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

namespace CppAD {

    template <class Base>
    inline CG<Base> pow(const CG<Base> &x, const CG<Base> &y) {
        if (x.isParameter() && y.isParameter()) {
            return CG<Base> (pow(x.getValue(), y.getValue()));
        }

        CodeHandler<Base>* handler;
        if (y.isParameter()) {
            if (y.IdenticalZero()) {
                return CG<Base> (Base(1.0)); // does not consider that x could be infinity
            } else if (y.IdenticalOne()) {
                return CG<Base> (x);
            }
            handler = x.getCodeHandler();
        } else {
            handler = y.getCodeHandler();
        }

        CG<Base> result(*handler, new SourceCodeFragment<Base>(CGPowOp, x.argument(), y.argument()));
        if (x.isValueDefined() && y.isValueDefined()) {
            result.setValue(pow(x.getValue(), y.getValue()));
        }
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
        if (x.isParameter()) {
            if (x.getValue() > Base(0.0)) {
                return CG<Base> (Base(1.0));
            } else if (x.getValue() == Base(0.0)) {
                return CG<Base> (Base(0.0));
            } else {
                return CG<Base> (Base(-1.0));
            }
        }

        CG<Base> result(*x.getCodeHandler(), new SourceCodeFragment<Base>(CGSignOp, x.argument()));
        if (x.isValueDefined()) {
            if (x.getValue() > Base(0.0)) {
                result.setValue(Base(1.0));
            } else if (x.getValue() == Base(0.0)) {
                result.setValue(Base(0.0));
            } else {
                result.setValue(Base(-1.0));
            }
        }
        return result;
    }

}

#endif

