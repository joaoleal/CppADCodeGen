/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppadcg/cg.hpp>

#include "test_solve.hpp"

bool SolveMul() {
    using namespace CppAD;
    using namespace std;

    typedef CG<double> CGD;

    // independent variable vector
    std::vector<AD<CGD> > u(2);
    u[0] = 4.0;
    u[1] = 2.0;
    
    Independent(u);

    // dependent variable vector
    std::vector< AD<CGD> > Z(4);

    // model
    Z[0] = u[0] * u[1];
    Z[1] = Z[0] * 4.; 
    Z[2] = 2. * Z[1]; 
    Z[3] = Z[2] * 1 - 64.0;
    
    // create f: U -> Z
    ADFun<CGD> fun(u, Z);

    bool ok = test_solve(fun, 3, 0, u);
    ok &= test_solve(fun, 3, 1, u);
    
    return  ok;
}

