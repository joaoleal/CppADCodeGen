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
