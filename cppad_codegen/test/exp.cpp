/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad_codegen/cppad_codegen.hpp>
#include <cmath>

#include "gcc_load_dynamic.hpp"

namespace { // BEGIN empty namespace

    bool ExpTestOne(void) {
        using CppAD::exp;
        using namespace CppAD;

        bool ok = true;

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(1);
        size_t s = 0;
        u[s] = 1.;

        CPPAD_TEST_VECTOR< AD<double> > U(1);
        U[s] = u[s];
        Independent(U);

        // dependent variable vector, indices, and values
        CPPAD_TEST_VECTOR< AD<double> > Z(2);
        size_t x = 0;
        size_t y = 1;
        Z[x] = exp(U[s]);
        Z[y] = exp(Z[x]);

        // define f : U -> Z and vectors for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("ExpTestOne", f, u, Z);

        return ok;
    }

    bool ExpTestTwo(void) {
        bool ok = true;

        using CppAD::exp;
        using namespace CppAD;

        // independent variable vector
        std::vector<double> u(1);
        size_t s = 0;
        u[s] = 1.;

        CPPAD_TEST_VECTOR< AD<double> > U(1);
        U[0] = u[s];
        Independent(U);

        // dependent variable vector 
        CPPAD_TEST_VECTOR< AD<double> > Z(1);
        Z[0] = exp(U[0]);

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("ExpTestTwo", f, u, Z);

        return ok;
    }

} // END empty namespace

bool Exp(void) {
    bool ok = true;
    ok &= ExpTestOne();
    ok &= ExpTestTwo();
    return ok;
}
