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

    bool LogTestOne(void) {
        using CppAD::log;
        using namespace CppAD;

        bool ok = true;

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(1);
        size_t s = 0;
        u[s] = 2.;

        CPPAD_TEST_VECTOR< AD<double> > U(1);
        U[s] = u[s];
        Independent(U);

        // dependent variable vector, indices, and values
        CPPAD_TEST_VECTOR< AD<double> > Z(2);
        size_t x = 0;
        size_t y = 1;
        Z[x] = log(U[s]);
        Z[y] = log(Z[x]);

        // define f : U -> Z and vectors for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("LogTestOne", f, u, Z, 1e-10, 1e-10);

        return ok;
    }

    bool LogTestTwo(void) {
        using CppAD::log;
        using namespace CppAD;

        bool ok = true;

        // independent variable vector
        std::vector<double> u(1);
        size_t s = 0;
        u[0] = 1.;

        CPPAD_TEST_VECTOR< AD<double> > U(1);
        U[0] = u[0];
        Independent(U);

        // a temporary values
        AD<double> x = exp(U[0]);

        // dependent variable vector 
        CPPAD_TEST_VECTOR< AD<double> > Z(1);
        Z[0] = log(x); // log( exp(u) )

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("LogTestTwo", f, u, Z, 1e-10, 1e-10);

        return ok;
    }

} // END empty namespace

bool Log(void) {
    bool ok = true;
    ok &= LogTestOne();
    ok &= LogTestTwo();
    return ok;
}
