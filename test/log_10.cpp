/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppadcg/cg.hpp>

#include "gcc_load_dynamic.hpp"
#include "log_10.hpp"

using namespace CppAD;

bool Log10() {
    // independent variable vector, indices, values, and declaration
    std::vector<double> u(1);
    size_t s = 0;
    u[s] = 10.;

    return test0nJac("log10", &Log10Func<double >, &Log10Func<CG<double> >, u, 1e-10, 1e-10);
}
