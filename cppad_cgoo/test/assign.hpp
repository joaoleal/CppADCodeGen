#ifndef CPPADCGOO_TEST_ASSIGN_INCLUDED
#define	CPPADCGOO_TEST_ASSIGN_INCLUDED
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
CppAD::ADFun<T>* AssignFunc(const std::vector<CppAD::AD<T> >& u) {
    using namespace CppAD;
    using namespace std;

    assert(u.size() == 2);

    std::vector<AD<T> > w(2);
    AD<T> a = u[0] + 2.0;
    AD<T> b = u[1] + 1.0;
    b = a;
    a += 5.0;
    w[0] = a;
    w[1] = b;

    // f(v) = |w|
    return new CppAD::ADFun<T > (u, w);
}

#endif