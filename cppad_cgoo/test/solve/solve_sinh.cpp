/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
#include <cmath>

#include <cppad_cgoo/cg.hpp>

#include "test_solve.hpp"

bool SolveSinh() {
    using namespace CppAD;
    using namespace std;

    typedef CG<double> CGD;

    // independent variable vector
    std::vector<AD<CGD> > u(2);
    u[0] = 2.0;
    u[1] = 1.0;

    Independent(u);

    // dependent variable vector
    std::vector< AD<CGD> > Z(2);

    // dependent variables
    Z[0] = sinh(u[0]);
    Z[1] = sinh(Z[0]) - u[1] * CppAD::sinh(CppAD::sinh(2.0));

    // create f: U -> Z
    ADFun<CGD> fun(u, Z);

    bool ok = test_solve(fun, 1, 0, u);

    return ok;
}
