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
#include "atan.hpp"

namespace { // BEGIN empty namespace

    bool AtanTestOne() {
        using namespace CppAD;

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(1);
        size_t s = 0;
        u[s] = 1.;

        bool ok = test0nJac("AtanTestOne", &AtanTestOneFunc<double >, &AtanTestOneFunc<CG<double> >, u);

        return ok;
    }

    bool AtanTestTwo() {
        using namespace CppAD;

        // independent variable vector
        std::vector<double> u(1);
        u[0] = 1.;

        bool ok = test0nJac("AtanTestTwo", &AtanTestTwoFunc<double >, &AtanTestTwoFunc<CG<double> >, u);

        return ok;
    }

} // END empty namespace

bool Atan(void) {
    bool ok = true;
    ok &= AtanTestOne();
    ok &= AtanTestTwo();
    return ok;
}
