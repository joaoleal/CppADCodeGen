#ifndef CPPADCGOO_TEST_MUL_INCLUDED
#define	CPPADCGOO_TEST_MUL_INCLUDED
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
CppAD::ADFun<T>* MulTestOneFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    // independent variable vector, indices, values, and declaration
    assert(U.size() == 2);

    size_t s = 0;
    size_t t = 1;

    // assign some parameters
    AD<T> zero = 0.;
    AD<T> one = 1.;

    // dependent variable vector and indices
    std::vector< AD<T> > Z(5);
    size_t x = 0;
    size_t y = 1;
    size_t z = 2;
    size_t u = 3;
    size_t v = 4;

    // assign the dependent variables
    Z[x] = U[s] * U[t]; // AD<double> * AD<double>
    Z[y] = Z[x] * 4.; // AD<double> *    double
    Z[z] = 4. * Z[y]; //    double  * AD<double> 
    Z[u] = one * Z[z]; // multiplication by parameter equal to one
    Z[v] = zero * Z[z]; // multiplication by parameter equal to zero

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

template<class T>
CppAD::ADFun<T>* MulTestTwoFunc(const std::vector<CppAD::AD<T> >& U) {
    using namespace CppAD;

    // independent variable vector
    assert(U.size() == 1);

    AD<T> a = U[0] * 1.; // AD<double> * double
    AD<T> b = a * 2; // AD<double> * int
    AD<T> c = 3. * b; // double     * AD<double> 
    AD<T> d = 4 * c; // int        * AD<double> 

    // dependent variable vector 
    std::vector< AD<T> > Z(1);
    Z[0] = U[0] * d; // AD<double> * AD<double>

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

#endif