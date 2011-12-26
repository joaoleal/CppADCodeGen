/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad/cppad_code_gen.hpp>

#include "gcc_load_dynamic.hpp"

namespace { // BEGIN empty namespace

    bool PowTestOne(void) {

        using CppAD::AD;
        using CppAD::NearEqual;

        bool ok = true;

        // domain space vector
        size_t n = 2;
        double x = 0.5;
        double y = 2.;
        std::vector<double> u(n);
        u[0] = x;
        u[1] = y;
        CPPAD_TEST_VECTOR< AD<double> > XY(n);
        XY[0] = u[0];
        XY[1] = u[1];

        // declare independent variables and start tape recording
        CppAD::Independent(XY);

        // range space vector 
        size_t m = 3;
        CPPAD_TEST_VECTOR< AD<double> > Z(m);
        Z[0] = CppAD::pow(XY[0], XY[1]); // pow(variable, variable)
        Z[1] = CppAD::pow(XY[0], y); // pow(variable, parameter)
        Z[2] = CppAD::pow(x, XY[1]); // pow(parameter, variable)

        // create f: XY -> Z and stop tape recording
        CppAD::ADFunCodeGen<double> f(XY, Z);

        ok &= test0nJac("PowTestOne", f, u, Z);

        return ok;
    }

    bool PowTestTwo(void) {
        bool ok = true;

        using CppAD::pow;
        using CppAD::exp;
        using namespace CppAD;


        // independent variable vector, indices, values, and declaration
        std::vector<double> uu(2);
        size_t s = 0;
        size_t t = 1;
        uu[s] = 2.;
        uu[t] = 3.;

        CPPAD_TEST_VECTOR< AD<double> > U(2);
        U[s] = uu[s];
        U[t] = uu[t];
        Independent(U);

        // dependent variable vector and indices
        CPPAD_TEST_VECTOR< AD<double> > Z(2);
        size_t x = 0;
        size_t y = 1;


        // dependent variable values
        AD<double> u = exp(U[s]); // u = exp(s)
        Z[x] = pow(u, U[t]); // x = exp(s * t)
        Z[y] = pow(Z[x], u); // y = exp( s * t * exp(s) )

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("PowTestTwo", f, uu, Z);

        return ok;
    }

    bool PowTestThree(void) {
        bool ok = true;

        using CppAD::AD;
        using CppAD::NearEqual;

        // domain space vector
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = 2.;
        CPPAD_TEST_VECTOR< AD<double> > x(n);
        x[0] = u[0];

        // declare independent variables and start tape recording
        CppAD::Independent(x);

        // range space vector 
        size_t m = 4;
        CPPAD_TEST_VECTOR< AD<double> > y(m);

        // some special cases
        y[0] = pow(x[0], 0.);
        y[1] = pow(0., x[0]);
        y[2] = pow(x[0], 1.);
        y[3] = pow(1., x[0]);

        // create f: x -> y and stop tape recording
        CppAD::ADFunCodeGen<double> f(x, y);

        ok &= test0nJac("PowTestThree", f, u, y);

        return ok;
    }

    bool PowTestFour(void) {
        bool ok = true;

        using CppAD::AD;
        using CppAD::NearEqual;

        // domain space vector
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = -2;
        CPPAD_TEST_VECTOR< AD<double> > x(n);
        x[0] = u[0];

        // declare independent variables and start tape recording
        CppAD::Independent(x);

        // range space vector 
        size_t m = 5;
        CPPAD_TEST_VECTOR< AD<double> > y(m);

        // some special cases (skip zero raised to a negative power)
        y[0] = pow(1., x[0]);
        size_t i;
        for (i = 1; i < m; i++)
            y[i] = pow(x[0], i - 1); // pow(AD<double>, int)

        // create f: x -> y and stop tape recording
        CppAD::ADFunCodeGen<double> f(x, y);

        ok &= test0nJac("PowTestFour", f, u, y);

        return ok;
    }

    bool PowTestFive(void) {
        bool ok = true;

        using CppAD::AD;
        using CppAD::NearEqual;

        // domain space vector
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = -1.;
        CPPAD_TEST_VECTOR< AD<double> > x(n);
        x[0] = u[0];

        // declare independent variables and start tape recording
        CppAD::Independent(x);

        // range space vector 
        size_t m = 1;
        CPPAD_TEST_VECTOR< AD<double> > y(m);

        // case of zero raised to a positive integer power
        double e = 2.;
        y[0] = pow(x[0], int(e)); // use pow(AD<double>, int)

        // create f: x -> y and stop tape recording
        CppAD::ADFunCodeGen<double> f(x, y);

        ok &= test0nJac("PowTestFive", f, u, y);

        return ok;
    }

    bool PowTestSix(void) {
        bool ok = true;

        using CppAD::AD;
        using CppAD::NearEqual;

        // domain space vector
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = 1.5;
        CPPAD_TEST_VECTOR< AD<double> > x(n);
        x[0] = u[0];

        // domain space vector
        CPPAD_TEST_VECTOR< AD<double> > X(n);
        X[0] = x[0];

        // declare independent variables and start tape recording
        CppAD::Independent(X);

        // range space vector 
        size_t m = 1;
        CPPAD_TEST_VECTOR< AD<double> > Y(m);

        // case of AD< AD<double> > raised to a double power
        double e = 2.5;
        Y[0] = pow(X[0], e);

        // create F: X -> Y and stop tape recording
        CppAD::ADFunCodeGen<double> f(X, Y);

        ok &= test0nJac("PowTestSix", f, u, Y);

        return ok;
    }

} // END empty namespace

bool Pow(void) {
    bool ok = true;
    ok &= PowTestOne();
    ok &= PowTestTwo();
    ok &= PowTestThree();
    ok &= PowTestFour();
    ok &= PowTestFive();
    ok &= PowTestSix();
    return ok;
}
