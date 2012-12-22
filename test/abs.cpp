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
#include "abs.hpp"

bool Abs() {
    using namespace CppAD;
    using namespace std;

    std::vector<std::vector<double> > uV;
    std::vector<double> u(1);
    u[0] = 0;
    uV.push_back(u);
    u[0] = 1;
    uV.push_back(u);
    u[0] = -1;
    uV.push_back(u);

    return test0nJac("abs", &AbsFunc<double >, &AbsFunc<CG<double> >, uV);
}
