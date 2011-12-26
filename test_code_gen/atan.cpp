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

    bool AtanTestOne() {
        using CppAD::atan;
        using namespace CppAD;

        bool ok = true;

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(1);
        size_t s = 0;
        u[s] = 1.;

        CPPAD_TEST_VECTOR< AD<double> > U(1);
        U[s] = u[s];
        Independent(U);

        // some temporary values
        AD<double> x = cos(U[s]);
        AD<double> y = sin(U[s]);
        AD<double> z = y / x; // tan(s)

        // dependent variable vector and indices
        CPPAD_TEST_VECTOR< AD<double> > Z(1);
        size_t a = 0;

        // dependent variable values
        Z[a] = atan(z); // atan( tan(s) )

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("AtanTestOne", f, u, Z);

        return ok;
    }

    bool AtanTestTwo(void) {
        using CppAD::atan;
        using CppAD::sin;
        using CppAD::cos;
        using namespace CppAD;

        bool ok = true;

        // independent variable vector
        std::vector<double> u(1);
        u[0] = 1.;
        CPPAD_TEST_VECTOR< AD<double> > U(1);
        U[0] = u[0];
        Independent(U);

        // a temporary values
        AD<double> x = sin(U[0]) / cos(U[0]);

        // dependent variable vector 
        CPPAD_TEST_VECTOR< AD<double> > Z(1);
        Z[0] = atan(x); // atan( tan(u) )

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("AtanTestTwo", f, u, Z);

        return ok;
    }

} // END empty namespace

bool Atan(void) {
    bool ok = true;
    ok &= AtanTestOne();
    ok &= AtanTestTwo();
    return ok;
}
