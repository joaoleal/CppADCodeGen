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
#include "div.hpp"

namespace { // BEGIN empty namespace

    bool DivTestOne() {
        bool ok = true;

        using namespace CppAD;

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(2);
        size_t s = 0;
        size_t t = 1;
        u[s] = 2.;
        u[t] = 3.;

        ok &= test0nJac("DivTestOne", &DivTestOneFunc<double >, &DivTestOneFunc<CG<double> >, u);

        return ok;
    }

    bool DivTestTwo() {
        using namespace CppAD;

        bool ok = true;

        // independent variable vector
        std::vector<double> u(1);
        u[0] = .5;

        ok &= test0nJac("DivTestTwo", &DivTestTwoFunc<double >, &DivTestTwoFunc<CG<double> >, u);

        return ok;
    }

    bool DivTestThree() {
        using namespace CppAD;

        bool ok = true;

        // more testing of variable / variable case 
        std::vector<double> u(2);
        u[0] = 2.;
        u[1] = 3.;

        ok &= test0nJac("DivTestThree", &DivTestThreeFunc<double >, &DivTestThreeFunc<CG<double> >, u);

        return ok;
    }

} // END empty namespace

bool Div(void) {
    bool ok = true;
    ok &= DivTestOne();
    ok &= DivTestTwo();
    ok &= DivTestThree();
    return ok;
}
