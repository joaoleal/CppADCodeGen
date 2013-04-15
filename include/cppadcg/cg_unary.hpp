#ifndef CPPAD_CG_UNARY_INCLUDED
#define CPPAD_CG_UNARY_INCLUDED
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
    inline CG<Base> CG<Base>::operator+() const {
        return CG<Base> (*this); // nothing to do
    }

    template<class Base>
    inline CG<Base> CG<Base>::operator-() const {
        if (isParameter()) {
            return CG<Base> (-getValue());

        } else {
            CG<Base> result(*getCodeHandler(), new SourceCodeFragment<Base>(CGUnMinusOp, this->argument()));
            if (isValueDefined()) {
                result.setValue(-getValue());
            }
            return result;
        }
    }

}

#endif

