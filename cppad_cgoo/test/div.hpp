#ifndef CPPADCGOO_TEST_DIV_INCLUDED
#define	CPPADCGOO_TEST_DIV_INCLUDED
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
CppAD::ADFun<T>* DivTestOneFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    assert(U.size() == 2);
    size_t s = 0;
    size_t t = 1;

    // assign some parameters
    AD<T> zero = 0.;
    AD<T> one = 1.;

    // dependent variable vector and indices
    std::vector< AD<T> > Z(6);
    size_t x = 0;
    size_t y = 1;
    size_t z = 2;
    size_t u = 3;
    size_t v = 4;
    size_t w = 5;

    // dependent variables
    Z[x] = U[s] / U[t]; // AD<double> / AD<double>
    Z[y] = Z[x] / 4.; // AD<double> / double
    Z[z] = 5. / Z[y]; //     double / AD<double> 
    Z[u] = Z[z] / one; // division by a parameter equal to one
    Z[v] = Z[z] / 1.; // division by a double equal to one
    Z[w] = zero / Z[z]; // division into a parameter equal to zero

    // create f : U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

template<class T>
CppAD::ADFun<T>* DivTestTwoFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    assert(U.size() == 1);

    // independent variable vector
    AD<T> a = U[0] / 1.; // AD<double> / double
    AD<T> b = a / 2; // AD<double> / int
    AD<T> c = 3. / b; // double     / AD<double> 
    AD<T> d = 4 / c; // int        / AD<double> 

    // dependent variable vector 
    std::vector< AD<T> > Z(1);
    Z[0] = U[0] * U[0] / d; // AD<double> / AD<double>

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

template<class T>
CppAD::ADFun<T>* DivTestThreeFunc(const std::vector<CppAD::AD<T> >& X) {
    using namespace CppAD;

    assert(X.size() == 2);

    size_t m = 1;
    std::vector< AD<T> > Y(m);
    Y[0] = X[0] / X[1];
    return new ADFun<T > (X, Y);
}

#endif