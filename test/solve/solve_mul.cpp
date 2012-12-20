/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

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

