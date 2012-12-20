/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

#include <cppadcg/cg.hpp>

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
