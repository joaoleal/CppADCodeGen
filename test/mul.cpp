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
