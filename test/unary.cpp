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
#include "unary.hpp"

using namespace CppAD;

namespace { // BEGIN empty namespace

    bool Plus() {
        std::vector<double> u(2);
        u[0] = -3.;
        u[1] = 2.;

        return test0nJac("UnaryPlus", &UnaryPlusFunc<double >, &UnaryPlusFunc<CG<double> >, u);
    }

    bool Minus() {
        std::vector<double> u(2);
        u[0] = -3.;
        u[1] = 2.;

        return test0nJac("UnaryMinus", &UnaryMinusFunc<double >, &UnaryMinusFunc<CG<double> >, u);
    }

} // END empty namespace

bool Unary() {
    bool ok = true;
    ok &= Plus();
    ok &= Minus();
    return ok;
}
