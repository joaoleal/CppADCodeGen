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
