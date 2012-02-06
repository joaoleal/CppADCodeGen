/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad_cgoo/cg.hpp>

#include "gcc_load_dynamic.hpp"
#include "assign.hpp"

bool Assign() {
    using namespace CppAD;
    using namespace std;

    std::vector<double> u(2);
    u[0] = 0;
    u[1] = 1;
    
    return test0nJac("assign", &AssignFunc<double >, &AssignFunc<CG<double> >, u);
}
