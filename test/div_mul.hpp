#ifndef CPPADCGOO_TEST_DIV_MUL_INCLUDED
#define	CPPADCGOO_TEST_DIV_MUL_INCLUDED
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

#include <assert.h>

template<class T>
CppAD::ADFun<T>* DivMulTestOneFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    assert(U.size() == 3);

    // dependent variable vector and indices
    std::vector< AD<T> > Z(1);

    // dependent variables
    Z[0] = U[0] / (U[1] * U[2]); // AD<double> / (AD<double> * AD<double>)

    // create f : U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

template<class T>
CppAD::ADFun<T>* DivMulTestTwoFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    assert(U.size() == 4);

    // dependent variable vector 
    std::vector< AD<T> > Z(1);
    Z[0] = U[0] / (U[1] * U[2]) / U[3]; // AD<double> / (AD<double> * AD<double>)

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

#endif