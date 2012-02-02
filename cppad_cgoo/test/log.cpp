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
#include "log.hpp"

using namespace CppAD;

namespace { // BEGIN empty namespace

    bool LogTestOne() {
        bool ok = true;

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(1);
        u[0] = 2.;

        ok &= test0nJac("LogTestOne", &LogTestOneFunc<double >, &LogTestOneFunc<CG<double> >, u, 1e-10, 1e-10);
        return ok;
    }

    bool LogTestTwo(void) {
        bool ok = true;

        // independent variable vector
        std::vector<double> u(1);
        u[0] = 1.;

        ok &= test0nJac("LogTestTwo", &LogTestTwoFunc<double >, &LogTestTwoFunc<CG<double> >, u, 1e-10, 1e-10);

        return ok;
    }

} // END empty namespace

bool Log() {
    bool ok = true;
    ok &= LogTestOne();
    ok &= LogTestTwo();
    return ok;
}
