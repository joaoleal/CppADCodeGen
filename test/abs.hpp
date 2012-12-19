#ifndef CPPADCGOO_TEST_ABS_INCLUDED
#define	CPPADCGOO_TEST_ABS_INCLUDED
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
CppAD::ADFun<T>* AbsFunc(const std::vector<CppAD::AD<T> >& u) {
    using namespace CppAD;
    using namespace std;

    assert(u.size() == 1);

    std::vector<CppAD::AD<T> > w(1);
    w[0] = CppAD::abs(u[0]);

    // f(v) = |w|
    return new CppAD::ADFun<T>(u, w);
}

#endif