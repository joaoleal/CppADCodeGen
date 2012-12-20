#ifndef CPPAD_CG_ORDERED_INCLUDED
#define CPPAD_CG_ORDERED_INCLUDED
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

