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
