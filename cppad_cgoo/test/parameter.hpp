#ifndef CPPADCGOO_TEST_PARAMETER_INCLUDED
#define	CPPADCGOO_TEST_PARAMETER_INCLUDED
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
CppAD::ADFun<T>* ParameterFunc(const std::vector<CppAD::AD<T> >& ax) {
    using namespace CppAD;

    // number of different parameter values
    size_t n_parameter = 7;

    // number of parameter repeats
    size_t n_repeat = 5;

    // independent variable vector
    size_t j, n = ax.size();
    assert(n == n_parameter * n_repeat);

    // dependent variable vector and indices
    size_t i, m = n;
    std::vector< AD<T> > ay(m);
    for (i = 0; i < m; i++) { // must avoid Float(k) = 0 because it would get optimized out	
        size_t k = (i % n_parameter);
        k = k * k * 10 + 1;
        j = i;
        ay[i] = ax[j] + T(k);
    }

    // create f: ax -> ay 
    return new ADFun<T > (ax, ay);
}

#endif