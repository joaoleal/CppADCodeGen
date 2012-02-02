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
#include "sub.hpp"

using namespace CppAD;

namespace { // BEGIN empty namespace

    bool One() {
        std::vector<double> u(2);
        size_t s = 0;
        size_t t = 1;
        u[s] = 3.;
        u[t] = 2.;

        return test0nJac("SubOne", &OneFunc<double >, &OneFunc<CG<double> >, u);
    }

    bool Two() {
        std::vector<double> u(1);
        u[0] = .5;

        return test0nJac("SubTwo", &TwoFunc<double >, &TwoFunc<CG<double> >, u, 1e-10, 1e-10);
    }

    bool Three() {
        std::vector<double> u(1);
        u[0] = 1.;

        return test0nJac("SubThree", &ThreeFunc<double >, &ThreeFunc<CG<double> >, u, 1e-10, 1e-10);
    }

    bool Four() {
        std::vector<double> u(1);
        u[0] = 1.;

        return test0nJac("SubFour", &FourFunc<double >, &FourFunc<CG<double> >, u, 1e-10, 1e-10);
    }


} // END empty namespace

bool Sub() {
    bool ok = true;
    ok &= One();
    ok &= Two();
    ok &= Three();
    ok &= Four();
    return ok;
}
