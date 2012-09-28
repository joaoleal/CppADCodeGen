#ifndef CPPADCGOO_TEST_DIV_MUL_INCLUDED
#define	CPPADCGOO_TEST_DIV_MUL_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

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