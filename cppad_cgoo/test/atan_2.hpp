#ifndef CPPADCGOO_TEST_ATAN_2_INCLUDED
#define	CPPADCGOO_TEST_ATAN_2_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <assert.h>

template<class T>
CppAD::ADFun<T>* Atan2Func(const std::vector<CppAD::AD<T> >& u) {
    using CppAD::atan;
    using CppAD::sin;
    using CppAD::cos;
    using namespace CppAD;

    assert(u.size() == 1);

    // a temporary values
    AD<T> x = cos(u[0]);
    AD<T> y = sin(u[0]);

    // dependent variable vector 
    std::vector< AD<T> > Z(1);
    Z[0] = atan2(y, x);

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (u, Z);
}

#endif