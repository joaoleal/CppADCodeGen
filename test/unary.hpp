#ifndef CPPADCGOO_TEST_UNARY_INCLUDED
#define	CPPADCGOO_TEST_UNARY_INCLUDED
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
CppAD::ADFun<T>* UnaryPlusFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    // independent variable vector, indices, values, and declaration
    assert(U.size() == 2);

    // dependent variable vector and indices
    std::vector< AD<T> > Z(2);

    // dependent variable values
    Z[0] = +U[0]; // + AD<double>
    Z[1] = +U[1]; // + AD<double>

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

template<class T>
CppAD::ADFun<T>* UnaryMinusFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    // independent variable vector, indices, values, and declaration
    assert(U.size() == 2);

    // dependent variable vector and indices
    std::vector< AD<T> > Z(2);

    // dependent variable values
    Z[0] = -U[0]; // + AD<double>
    Z[1] = -U[1]; // + AD<double>

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}


#endif
