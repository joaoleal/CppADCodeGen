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

bool Atan2() {
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
    AD<double> x = cos(U[0]);
    AD<double> y = sin(U[0]);

    // dependent variable vector 
    CPPAD_TEST_VECTOR< AD<double> > Z(1);
    Z[0] = atan2(y, x);

    // create f: U -> Z and vectors used for derivative calculations
    ADFunCodeGen<double> f(U, Z);

    ok &= test0nJac("Atan2", f, u, Z);

    return ok;
}
