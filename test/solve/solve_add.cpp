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

bool SolveAdd() {
    using namespace CppAD;
    using namespace std;

    typedef CG<double> CGD;

    // independent variable vector
    std::vector<AD<CGD> > u(2);
    u[0] = -1.5;
    u[1] = -0.5;
    Independent(u);

    // dependent variable vector
    std::vector< AD<CGD> > Z(3);

    // model
    Z[0] = u[0] + u[1]; // AD<double> + AD<double>
    Z[1] = Z[0] + 1.; // AD<double> + double
    Z[2] = 1. + Z[1]; // double + AD<double> 

    // create f: U -> Z 
    ADFun<CGD> fun(u, Z);

    bool ok = test_solve(fun, 2, 1, u);
    
    return ok;
}
