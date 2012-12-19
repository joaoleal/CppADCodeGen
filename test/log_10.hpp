#ifndef CPPADCGOO_TEST_LOG_10_INCLUDED
#define	CPPADCGOO_TEST_LOG_10_INCLUDED
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
CppAD::ADFun<T>* Log10Func(const std::vector<CppAD::AD<T> >& U) {
    using CppAD::log10;
    using CppAD::log;
    using namespace CppAD;

    assert(U.size() == 1);

    // dependent variable vector, indices, and values
    std::vector< AD<T> > Z(2);
    size_t x = 0;
    size_t y = 1;
    Z[x] = log10(U[0]);
    Z[y] = log10(Z[x]);

    // define f : U -> Z and vectors for derivative calculations
    return new ADFun<T > (U, Z);
}

#endif