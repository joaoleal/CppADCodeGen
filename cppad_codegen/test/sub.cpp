/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad_codegen/cppad_codegen.hpp>

#include "gcc_load_dynamic.hpp"

namespace { // BEGIN empty namespace

    bool One() {
        using namespace CppAD;

        bool ok = true;

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(2);
        size_t s = 0;
        size_t t = 1;
        u[s] = 3.;
        u[t] = 2.;

        CPPAD_TEST_VECTOR< AD<double> > U(2);
        U[s] = u[0];
        U[t] = u[1];
        Independent(U);

        // dependent variable vector and indices
        CPPAD_TEST_VECTOR< AD<double> > Z(3);
        size_t x = 0;
        size_t y = 1;
        size_t z = 2;

        // dependent variable values
        Z[x] = U[s] - U[t]; // AD<double> - AD<double>
        Z[y] = Z[x] - 1.; // AD<double> - double
        Z[z] = 1. - Z[y]; // double - AD<double> 

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("SubOne", f, u, Z);

        return ok;

    }

    bool Two() {
        using namespace CppAD;

        bool ok = true;

        // independent variable vector
        std::vector<double> u(1);
        u[0] = .5;
        CPPAD_TEST_VECTOR< AD<double> > U(1);
        U[0] = u[0];
        Independent(U);

        AD<double> a = 2. * U[0] - 1.; // AD<double> - double
        AD<double> b = a - 2; // AD<double> - int
        AD<double> c = 3. - b; // double     - AD<double> 
        AD<double> d = 4 - c; // int        - AD<double> 

        // dependent variable vector 
        CPPAD_TEST_VECTOR< AD<double> > Z(1);
        Z[0] = U[0] - d; // AD<double> - AD<double>

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("SubTwo", f, u, Z, 1e-10, 1e-10);

        return ok;
    }

    bool Three() {
        bool ok = true;
        using namespace CppAD;

        // special cases where tests above check OK and SubpvOp 
        // implementation is known to be worng. 
        // Probably two minuses make a plus.
        size_t n = 1;
        std::vector<double> u(1);
        u[0] = 1.;
        CPPAD_TEST_VECTOR< AD<double> > X(n);
        X[0] = u[0];
        Independent(X);
        size_t m = 1;
        CPPAD_TEST_VECTOR< AD<double> > Y(m);
        Y[0] = 1. - X[0];
        ADFunCodeGen<double> f(X, Y);

        ok &= test0nJac("SubThree", f, u, Y, 1e-10, 1e-10);

        return ok;
    }

    bool Four() {
        using namespace CppAD;

        bool ok = true;

        // special cases where parameter number is equal to
        // variable index in result.
        size_t n = 1;
        std::vector<double> u(1);
        u[0] = 1.;
        CPPAD_TEST_VECTOR< AD<double> > X(n);
        X[0] = u[0];
        Independent(X);
        size_t m = 1;
        CPPAD_TEST_VECTOR< AD<double> > Y(m);
        if (0. < X[0] && X[0] < 10.)
            Y[0] = X[0] - 2.;
        else Y[0] = X[0] - 2.;
        ADFunCodeGen<double> f(X, Y);

        ok &= test0nJac("SubFour", f, u, Y, 1e-10, 1e-10);

        return ok;
    }


} // END empty namespace

bool Sub() {
    bool ok = true;
    ok &= One();
    ok &= Two();
    ok &= Three();
    ok &= Four();
    return ok;
}
