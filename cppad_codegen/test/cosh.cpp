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

bool Cosh(void) {
    bool ok = true;

    using CppAD::sinh;
    using CppAD::cosh;
    using namespace CppAD;

    // independent variable vector
    std::vector<double> u(1);
    u[0] = 1.;
    CPPAD_TEST_VECTOR< AD<double> > U(1);
    U[0] = u[0];
    Independent(U);

    // dependent variable vector 
    CPPAD_TEST_VECTOR< AD<double> > Z(1);
    Z[0] = cosh(U[0]);

    // create f: U -> Z and vectors used for derivative calculations
    ADFunCodeGen<double> f(U, Z);

    ok &= test0nJac("cosh", f, u, Z, 1e-10, 1e-10);

    return ok;
}
