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
#include "exp.hpp"

using namespace CppAD;

namespace { // BEGIN empty namespace

    bool ExpTestOne() {
        bool ok = true;

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(1);
        size_t s = 0;
        u[s] = 1.;

        ok &= test0nJac("ExpTestOne", &ExpTestOneFunc<double >, &ExpTestOneFunc<CG<double> >, u);

        return ok;
    }

    bool ExpTestTwo(void) {
        bool ok = true;

        // independent variable vector
        std::vector<double> u(1);
        size_t s = 0;
        u[s] = 1.;

        ok &= test0nJac("ExpTestTwo", &ExpTestTwoFunc<double >, &ExpTestTwoFunc<CG<double> >, u);

        return ok;
    }

} // END empty namespace

bool Exp() {
    bool ok = true;
    ok &= ExpTestOne();
    ok &= ExpTestTwo();
    return ok;
}
