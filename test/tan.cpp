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
#include "tan.hpp"

using namespace CppAD;

namespace {

    bool tan() {
        double eps = 100. * std::numeric_limits<double>::epsilon();

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(1);
        u[0] = .7;

        bool ok = test0nJac("tan_first", &tanFirstFunc<double >, &tanFirstFunc<CG<double> >, u, eps, eps);
        
        ok &= test0nJac("tan_last", &tanLastFunc<double >, &tanLastFunc<CG<double> >, u, eps, eps);
        
        return ok;
    }

    bool tanh() {
        double eps = 100. * std::numeric_limits<double>::epsilon();

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(1);
        u[0] = .5;

        bool ok = test0nJac("tanh_first", &tanhFirstFunc<double >, &tanhFirstFunc<CG<double> >, u, eps, eps);
        
        ok &= test0nJac("tanh_last", &tanhLastFunc<double >, &tanhLastFunc<CG<double> >, u, eps, eps);
        
        return ok;
    }
}

bool Tan() {
    bool ok = true;
    ok &= tan();
    ok &= tanh();
    return ok;
}
