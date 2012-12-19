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
#include "parameter.hpp"

bool Parameter() {
    using namespace CppAD;

    // number of different parameter values
    size_t n_parameter = 7;

    // number of parameter repeats
    size_t n_repeat = 5;

    // independent variable vector
    size_t n = n_parameter * n_repeat;
    std::vector<double> u(n);
    for (size_t j = 0; j < n; j++) {
        u[j] = double(j);
    }

    return test0nJac("parameter", &ParameterFunc<double >, &ParameterFunc<CG<double> >, u);
}

