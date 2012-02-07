#ifndef CPPADCGOO_TEST_SUB_INCLUDED
#define	CPPADCGOO_TEST_SUB_INCLUDED
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
CppAD::ADFun<T>* OneFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    // independent variable vector, indices, values, and declaration
    assert(U.size() == 2);
    size_t s = 0;
    size_t t = 1;

    // dependent variable vector and indices
    std::vector< AD<T> > Z(3);
    size_t x = 0;
    size_t y = 1;
    size_t z = 2;

    // dependent variable values
    Z[x] = U[s] - U[t]; // AD<double> - AD<double>
    Z[y] = Z[x] - 1.; // AD<double> - double
    Z[z] = 1. - Z[y]; // double - AD<double> 

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

template<class T>
CppAD::ADFun<T>* TwoFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    assert(U.size() == 1);

    AD<T> a = 2. * U[0] - 1.; // AD<double> - double
    AD<T> b = a - 2; // AD<double> - int
    AD<T> c = 3. - b; // double     - AD<double> 
    AD<T> d = 4 - c; // int        - AD<double> 

    // dependent variable vector 
    std::vector< AD<T> > Z(1);
    Z[0] = U[0] - d; // AD<double> - AD<double>

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

template<class T>
CppAD::ADFun<T>* ThreeFunc(const std::vector<CppAD::AD<T> >& X) {
    using namespace CppAD;

    assert(X.size() == 1);
    // special cases where tests above check OK and SubpvOp 
    // implementation is known to be wrong. 
    // Probably two minuses make a plus.

    size_t m = 1;
    std::vector< AD<T> > Y(m);
    Y[0] = 1. - X[0];

    return new ADFun<T > (X, Y);
}

template<class T>
CppAD::ADFun<T>* FourFunc(const std::vector<CppAD::AD<T> >& X) {
    using namespace CppAD;

    assert(X.size() == 1);

    // special cases where parameter number is equal to
    // variable index in result.
    size_t m = 1;
    std::vector< AD<T> > Y(m);
    Y[0] = X[0] - 2.;

    return new ADFun<T > (X, Y);
}

#endif
