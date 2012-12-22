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
