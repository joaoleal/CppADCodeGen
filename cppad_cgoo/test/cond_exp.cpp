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
#include "cond_exp.hpp"

using namespace CppAD;
using namespace std;

namespace { // Begin empty namespace

    bool CondExp_pvvv() {
        bool ok = true;

        // independent variable vector
        std::vector<double> u(3);
        u[0] = 0.;
        u[1] = 1.;
        u[2] = 2.;

        ok &= test0nJac("CondExp_pvvv", &CondExp_pvvvFunc<double >, &CondExp_pvvvFunc<CG<double> >, u);

        return ok;
    }

    bool CondExp_vpvv() {
        bool ok = true;

        // independent variable vector
        std::vector<double> u(3);
        u[0] = 0.;
        u[1] = 1.;
        u[2] = 2.;

        ok &= test0nJac("CondExp_vpvv", &CondExp_vpvvFunc<double >, &CondExp_vpvvFunc<CG<double> >, u);

        return ok;
    }

    bool CondExp_vvpv() {
        bool ok = true;

        // independent variable vector
        std::vector<double> u(3);
        u[0] = 0.;
        u[1] = 1.;
        u[2] = 2.;

        ok &= test0nJac("CondExp_vvpv", &CondExp_vvpvFunc<double >, &CondExp_vvpvFunc<CG<double> >, u);

        return ok;
    }

    bool CondExp_vvvp() {
        bool ok = true;

        // independent variable vector
        std::vector<double> u(3);
        u[0] = 0.;
        u[1] = 1.;
        u[2] = 2.;

        ok &= test0nJac("CondExp_vvvp", &CondExp_vvvpFunc<double >, &CondExp_vvvpFunc<CG<double> >, u);

        return ok;
    }

} // end empty namespace

bool CondExp() {
    bool ok = true;
    ok &= CondExp_pvvv();
    ok &= CondExp_vpvv();
    ok &= CondExp_vvpv();
    ok &= CondExp_vvvp();
    return ok;
}
// END PROGRAM
