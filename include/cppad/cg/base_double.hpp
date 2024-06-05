#ifndef CPPAD_CG_BASE_DOUBLE_INCLUDED
#define CPPAD_CG_BASE_DOUBLE_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

namespace CppAD {

/**
 * Specialization of the numeric_limits for doubles
 */
CPPAD_NUMERIC_LIMITS(double, cg::CG<double>)

/**
 * Specialization of the machine epsilon for CG<double>
 */
template <>
inline cg::CG<double> epsilon<cg::CG<double> >() {
    return std::numeric_limits<double>::epsilon();
}

/**
 * Absolute Zero multiplication
 */
inline cg::CG<double> azmul(const cg::CG<double>& x, const cg::CG<double>& y) {
    cg::CG<double> zero(0.0);
    if (x == zero)
        return zero;
    return x * y;
}
//CPPAD_AZMUL( cg::CG<double> )

} // END CppAD namespace

#endif

