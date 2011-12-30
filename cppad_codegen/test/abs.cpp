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

bool Abs() {
    using namespace CppAD;
    using namespace std;

    bool ok = true;

    std::vector<double> u(1);
    u[0] = 0;

    CPPAD_TEST_VECTOR<CppAD::AD<double> > U(1);
    U[0] = u[0];
    Independent(U);

    CPPAD_TEST_VECTOR<CppAD::AD<double> > w(1);
    w[0] = CppAD::abs(U[0]);

    // f(v) = |w|
    CppAD::ADFunCodeGen<double> f(U, w);

    ok &= test0nJac("abs", f, u, w);

    return ok;
}
