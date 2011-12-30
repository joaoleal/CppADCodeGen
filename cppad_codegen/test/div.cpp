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

    bool DivTestOne(void) {
        bool ok = true;

        using namespace CppAD;

        // assign some parameters
        AD<double> zero = 0.;
        AD<double> one = 1.;

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
        CPPAD_TEST_VECTOR< AD<double> > Z(6);
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

        // check division into a zero valued parameter results in a parameter
        // (must do this before creating f because it erases the tape)
        ok &= Parameter(Z[w]);

        // create f : U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        // check parameter flag
        ok &= f.Parameter(w);

        ok &= test0nJac("DivTestOne", f, uu, Z);

        return ok;
    }

    bool DivTestTwo(void) {
        using namespace CppAD;

        bool ok = true;

        // independent variable vector
        std::vector<double> u(1);
        u[0] = .5;
        CPPAD_TEST_VECTOR< AD<double> > U(1);
        U[0] = u[0];
        Independent(U);

        AD<double> a = U[0] / 1.; // AD<double> / double
        AD<double> b = a / 2; // AD<double> / int
        AD<double> c = 3. / b; // double     / AD<double> 
        AD<double> d = 4 / c; // int        / AD<double> 

        // dependent variable vector 
        CPPAD_TEST_VECTOR< AD<double> > Z(1);
        Z[0] = U[0] * U[0] / d; // AD<double> / AD<double>

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("DivTestTwo", f, u, Z);

        return ok;
    }

    bool DivTestThree(void) {
        using namespace CppAD;

        bool ok = true;

        // more testing of variable / variable case 
        std::vector<double> u(2);
        u[0] = 2.;
        u[1] = 3.;
        size_t n = 2;
        CPPAD_TEST_VECTOR< AD<double> > X(n);
        X[0] = u[0];
        X[1] = u[1];
        Independent(X);
        size_t m = 1;
        CPPAD_TEST_VECTOR< AD<double> > Y(m);
        Y[0] = X[0] / X[1];
        ADFunCodeGen<double> f(X, Y);

        ok &= test0nJac("DivTestThree", f, u, X);

        return ok;
    }

} // END empty namespace

bool Div(void) {
    bool ok = true;
    ok &= DivTestOne();
    ok &= DivTestTwo();
    return ok;
}
