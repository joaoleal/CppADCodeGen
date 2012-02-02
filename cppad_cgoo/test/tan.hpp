#ifndef CPPADCGOO_TEST_TAN_INCLUDED
#define	CPPADCGOO_TEST_TAN_INCLUDED
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
CppAD::ADFun<T>* tanFirstFunc(const std::vector<CppAD::AD<T> >& ax) {
    using CppAD::AD;

    assert(ax.size() == 1);

    // dependent variable vector and indices
    std::vector< AD<T> > ay(1);

    ay[0] = atan(tan(ax[0]));

    // create f: x -> y and vectors used for derivative calculations
    return new CppAD::ADFun<T > (ax, ay);
}

template<class T>
CppAD::ADFun<T>* tanLastFunc(const std::vector<CppAD::AD<T> >& ax) {
    using CppAD::AD;

    assert(ax.size() == 1);

    // dependent variable vector and indices
    std::vector< AD<T> > ay(1);

    ay[0] = tan(atan(ax[0]));

    // create f: x -> y and vectors used for derivative calculations
    return new CppAD::ADFun<T > (ax, ay);
}

template<class T>
CppAD::ADFun<T>* tanhFirstFunc(const std::vector<CppAD::AD<T> >& ax) {
    using CppAD::AD;

    // independent variable vector, indices, values, and declaration
    assert(ax.size() == 1);

    // dependent variable vector and indices
    std::vector< AD<T> > ay(1);
    AD<T> z = tanh(ax[0]);
    ay[0] = .5 * log((1. + z) / (1. - z));

    // create f: x -> y and vectors used for derivative calculations
    return new CppAD::ADFun<T > (ax, ay);
}

template<class T>
CppAD::ADFun<T>* tanhLastFunc(const std::vector<CppAD::AD<T> >& ax) {
    using CppAD::AD;

    // independent variable vector, indices, values, and declaration
    assert(ax.size() == 1);

    // dependent variable vector and indices
    std::vector< AD<T> > ay(1);
    AD<T> z = .5 * log((1. + ax[0]) / (1. - ax[0]));
    ay[0] = tanh(z);

    // create f: x -> y and vectors used for derivative calculations
    return new CppAD::ADFun<T > (ax, ay);
}


#endif