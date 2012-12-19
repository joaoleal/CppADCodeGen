#ifndef CPPADCGOO_TEST_ATAN_INCLUDED
#define	CPPADCGOO_TEST_ATAN_INCLUDED
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
CppAD::ADFun<T>* AtanTestOneFunc(const std::vector<CppAD::AD<T> >& u) {
    using CppAD::atan;
    using namespace CppAD;

    assert(u.size() == 1);
    
    size_t s = 0;

    // some temporary values
    AD<T> x = cos(u[s]);
    AD<T> y = sin(u[s]);
    AD<T> z = y / x; // tan(s)

    // dependent variable vector and indices
    std::vector< AD<T> > Z(1);
    size_t a = 0;

    // dependent variable values
    Z[a] = atan(z); // atan( tan(s) )

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (u, Z);
}

template<class T>
CppAD::ADFun<T>* AtanTestTwoFunc(const std::vector<CppAD::AD<T> >& u) {
    using CppAD::atan;
    using CppAD::sin;
    using CppAD::cos;
    using namespace CppAD;

    assert(u.size() == 1);
    
    // a temporary values
    AD<T> x = sin(u[0]) / cos(u[0]);

    // dependent variable vector 
    std::vector< AD<T> > Z(1);
    Z[0] = atan(x); // atan( tan(u) )

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (u, Z);
}

#endif