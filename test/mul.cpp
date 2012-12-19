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
#include "mul.hpp"

using namespace CppAD;

namespace { // BEGIN empty namespace

    bool MulTestOne() {
        // independent variable vector, indices, values, and declaration
        std::vector<double> u(2);
        size_t s = 0;
        size_t t = 1;
        u[s] = 3.;
        u[t] = 2.;

        return test0nJac("MulTestOne", &MulTestOneFunc<double >, &MulTestOneFunc<CG<double> >, u);
    }

    bool MulTestTwo() {
        // independent variable vector
        std::vector<double> u(1);
        u[0] = .5;

        return test0nJac("MulTestTwo", &MulTestTwoFunc<double >, &MulTestTwoFunc<CG<double> >, u);
    }

} // END empty namespace

bool Mul() {
    bool ok = true;
    ok &= MulTestOne();
    ok &= MulTestTwo();
    return ok;
}
