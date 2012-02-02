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
#include "pow.hpp"

using namespace CppAD;

namespace { // BEGIN empty namespace

    bool PowTestOne() {
        // domain space vector
        size_t n = 2;
        double x = 0.5;
        double y = 2.;
        std::vector<double> u(n);
        u[0] = x;
        u[1] = y;

        return test0nJac("PowTestOne", &PowTestOneFunc<double >, &PowTestOneFunc<CG<double> >, u);
    }

    bool PowTestTwo() {
        // independent variable vector, indices, values, and declaration
        std::vector<double> u(2);
        size_t s = 0;
        size_t t = 1;
        u[s] = 2.;
        u[t] = 3.;

        return test0nJac("PowTestTwo", &PowTestTwoFunc<double >, &PowTestTwoFunc<CG<double> >, u);
    }

    bool PowTestThree() {
        // domain space vector
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = 2.;

        return test0nJac("PowTestThree", &PowTestThreeFunc<double >, &PowTestThreeFunc<CG<double> >, u);
    }

    bool PowTestFour() {
        // domain space vector
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = -2;

        return test0nJac("PowTestFour", &PowTestFourFunc<double >, &PowTestFourFunc<CG<double> >, u);
    }

    bool PowTestFive(void) {
        // domain space vector
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = -1.;

        return test0nJac("PowTestFive", &PowTestFiveFunc<double >, &PowTestFiveFunc<CG<double> >, u);
    }

    bool PowTestSix() {
        // domain space vector
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = 1.5;

        return test0nJac("PowTestSix", &PowTestSixFunc<double >, &PowTestSixFunc<CG<double> >, u);
    }

} // END empty namespace

bool Pow(void) {
    bool ok = true;
    ok &= PowTestOne();
    ok &= PowTestTwo();
    ok &= PowTestThree();
    ok &= PowTestFour();
    ok &= PowTestFive();
    ok &= PowTestSix();
    return ok;
}
