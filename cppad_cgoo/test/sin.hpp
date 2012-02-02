#ifndef CPPADCGOO_TEST_SIN_INCLUDED
#define	CPPADCGOO_TEST_SIN_INCLUDED
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
CppAD::ADFun<T>* SinFunc(const std::vector<CppAD::AD<T> >& U) {
    using CppAD::sin;
    using CppAD::cos;
    using namespace CppAD;

    // independent variable vector
    assert(U.size() == 1);

    // dependent variable vector 
    std::vector< AD<T> > Z(1);
    Z[0] = sin(U[0]);

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T>(U, Z);
}

#endif