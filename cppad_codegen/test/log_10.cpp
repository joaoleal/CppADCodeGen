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

bool Log10(void) {
    bool ok = true;

    using CppAD::log10;
    using CppAD::log;
    using namespace CppAD;

    // independent variable vector, indices, values, and declaration
    std::vector<double> u(1);
    size_t s = 0;
    u[s] = 10.;

    CPPAD_TEST_VECTOR< AD<double> > U(1);
    U[s] = u[s];
    Independent(U);

    // dependent variable vector, indices, and values
    CPPAD_TEST_VECTOR< AD<double> > Z(2);
    size_t x = 0;
    size_t y = 1;
    Z[x] = log10(U[s]);
    Z[y] = log10(Z[x]);

    // define f : U -> Z and vectors for derivative calculations
    ADFunCodeGen<double> f(U, Z);
  
    ok &= test0nJac("log10", f, u, Z, 1e-10, 1e-10);

    return ok;
}
