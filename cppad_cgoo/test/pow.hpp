#ifndef CPPADCGOO_TEST_POW_INCLUDED
#define	CPPADCGOO_TEST_POW_INCLUDED
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
CppAD::ADFun<T>* PowTestOneFunc(const std::vector<CppAD::AD<T> >& XY) {
    using CppAD::AD;

    assert(XY.size() == 2);

    double x = 0.5;
    double y = 2.;

    // range space vector 
    size_t m = 3;
    std::vector< AD<T> > Z(m);
    Z[0] = CppAD::pow(XY[0], XY[1]); // pow(variable, variable)
    Z[1] = CppAD::pow(XY[0], y); // pow(variable, parameter)
    Z[2] = CppAD::pow(x, XY[1]); // pow(parameter, variable)

    // create f: XY -> Z and stop tape recording
    return new CppAD::ADFun<T > (XY, Z);
}

template<class T>
CppAD::ADFun<T>* PowTestTwoFunc(const std::vector<CppAD::AD<T> >& U) {
    using CppAD::pow;
    using CppAD::exp;
    using namespace CppAD;

    assert(U.size() == 2);
    size_t s = 0;
    size_t t = 1;

    // dependent variable vector and indices
    std::vector< AD<T> > Z(2);
    size_t x = 0;
    size_t y = 1;

    // dependent variable values
    AD<T> u = exp(U[s]); // u = exp(s)
    Z[x] = pow(u, U[t]); // x = exp(s * t)
    Z[y] = pow(Z[x], u); // y = exp( s * t * exp(s) )

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<T > (U, Z);
}

template<class T>
CppAD::ADFun<T>* PowTestThreeFunc(const std::vector<CppAD::AD<T> >& x) {
    using CppAD::AD;

    // domain space vector
    assert(x.size() == 1);

    // range space vector 
    size_t m = 4;
    std::vector< AD<T> > y(m);

    // some special cases
    y[0] = pow(x[0], 0.);
    y[1] = pow(0., x[0]);
    y[2] = pow(x[0], 1.);
    y[3] = pow(1., x[0]);

    // create f: x -> y and stop tape recording
    return new CppAD::ADFun<T > (x, y);
}

template<class T>
CppAD::ADFun<T>* PowTestFourFunc(const std::vector<CppAD::AD<T> >& x) {
    using namespace CppAD;

    // domain space vector
    assert(x.size() == 1);

    // range space vector 
    int m = 5;
    std::vector< AD<T> > y(m);

    // some special cases (skip zero raised to a negative power)
    y[0] = pow(1., x[0]);
    for (int i = 1; i < m; i++)
        y[i] = pow(x[0], i - 1); // pow(AD<double>, int)

    // create f: x -> y and stop tape recording
    return new ADFun<T > (x, y);
}

template<class T>
CppAD::ADFun<T>* PowTestFiveFunc(const std::vector<CppAD::AD<T> >& x) {
    using CppAD::AD;

    // domain space vector
    assert(x.size() == 1);

    // range space vector 
    size_t m = 1;
    std::vector< AD<T> > y(m);

    // case of zero raised to a positive integer power
    double e = 2.;
    y[0] = pow(x[0], int(e)); // use pow(AD<double>, int)

    // create f: x -> y and stop tape recording
    return new CppAD::ADFun<T > (x, y);
}

template<class T>
CppAD::ADFun<T>* PowTestSixFunc(const std::vector<CppAD::AD<T> >& X) {
    using CppAD::AD;

    // domain space vector
    assert(X.size() == 1);

    // range space vector 
    size_t m = 1;
    std::vector< AD<T> > Y(m);

    // case of AD< AD<double> > raised to a double power
    double e = 2.5;
    Y[0] = pow(X[0], e);

    // create F: X -> Y and stop tape recording
    return new CppAD::ADFun<T>(X, Y);
}

#endif