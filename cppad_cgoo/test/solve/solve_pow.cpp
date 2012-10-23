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

bool SolvePow() {
    using namespace CppAD;
    using namespace std;

    typedef CG<double> CGD;

    // independent variable vector
    std::vector<AD<CGD> > u(4);
    u[0] = 4.0;
    u[1] = 2.0;
    u[2] = 1.0;
    u[3] = 4.0;

    Independent(u);

    // dependent variable vector
    std::vector< AD<CGD> > Z(3);

    // dependent variables
    Z[0] = pow(pow(u[0], u[1]), u[2]);
    Z[1] = Z[0] - u[2] * pow(4, 2);
    Z[2] = pow(Z[0], 2.0) - pow(pow(u[3], 2), 2);

    // create f: U -> Z
    ADFun<CGD> fun(u, Z);

    bool ok = true;//test_solve(fun, 1, 0, u);
    //ok &= test_solve(fun, 1, 1, u);
    ok &= test_solve(fun, 2, 0, u);

    return ok;
}
