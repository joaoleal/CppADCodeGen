#ifndef CPPADCGOO_TEST_ASIN_INCLUDED
#define	CPPADCGOO_TEST_ASIN_INCLUDED
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
CppAD::ADFun<T>* AsinFunc(const std::vector<CppAD::AD<T> >& u) {
    using namespace CppAD;
    using namespace std;
    using CppAD::asin;
    using namespace CppAD;

    // a temporary values
    AD<T> x = sin(u[0]);

    // dependent variable vector 
    std::vector< AD<T> > Z(1);
    Z[0] = asin(x); // asin( sin(u) )

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T>(u, Z);
}

#endif