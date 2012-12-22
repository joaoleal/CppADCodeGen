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
#include "cosh.hpp"

bool Cosh() {
    bool ok = true;

    using namespace CppAD;

    // independent variable vector
    std::vector<double> u(1);
    u[0] = 1.;

    ok &= test0nJac("cosh", &CoshFunc<double >, &CoshFunc<CG<double> >, u, 1e-10, 1e-10);

    return ok;
}
