#ifndef CPPADCGOO_TEST_EXP_INCLUDED
#define	CPPADCGOO_TEST_EXP_INCLUDED
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
CppAD::ADFun<T>* ExpTestOneFunc(const std::vector<CppAD::AD<T> >& U) {
    using CppAD::exp;
    using namespace CppAD;

    assert(U.size() == 1);
    size_t s = 0;

    // dependent variable vector, indices, and values
    std::vector< AD<T> > Z(2);
    size_t x = 0;
    size_t y = 1;
    Z[x] = exp(U[s]);
    Z[y] = exp(Z[x]);

    // define f : U -> Z and vectors for derivative calculations
    return new ADFun<T > (U, Z);
}

template<class T>
CppAD::ADFun<T>* ExpTestTwoFunc(const std::vector<CppAD::AD<T> >& U) {
    using CppAD::exp;
    using namespace CppAD;

    assert(U.size() == 1);

    // dependent variable vector 
    std::vector< AD<T> > Z(1);
    Z[0] = exp(U[0]);

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

#endif