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
#include "div_mul.hpp"

namespace { // BEGIN empty namespace

    bool DivMulTestOne() {
        bool ok = true;

        using namespace CppAD;

        // independent variable vector, indices, values, and declaration
        std::vector<double> u(3);
        u[0] = 2.;
        u[1] = 3.;
        u[2] = 4.;

        ok &= test0nJac("DivMulTestOne", &DivMulTestOneFunc<double >, &DivMulTestOneFunc<CG<double> >, u);

        return ok;
    }

    bool DivMulTestTwo() {
        using namespace CppAD;

        bool ok = true;

        // independent variable vector
        std::vector<double> u(4);
        u[0] = 2.;
        u[1] = 3.;
        u[2] = 4.;
        u[3] = 5.;

        ok &= test0nJac("DivMulTestTwo", &DivMulTestTwoFunc<double >, &DivMulTestTwoFunc<CG<double> >, u);

        return ok;
    }

} // END empty namespace

bool DivMul() {
    bool ok = true;
    ok &= DivMulTestOne();
    ok &= DivMulTestTwo();
    return ok;
}
