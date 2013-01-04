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

using namespace CppAD;

bool InputStream() {

    std::istringstream is("5.5"); // float
    CppAD::CG<double> var;

    is >> var;

    bool ok = var.isParameter();
    ok &= CppAD::NearEqual(var.getParameterValue(), 5.5, 10e-10, 10e-10);


    std::istringstream is2("8"); // integer
    CppAD::CG<double> var2;

    is2 >> var2;

    ok &= var2.isParameter();
    ok &= CppAD::NearEqual(var2.getParameterValue(), 8.0, 10e-10, 10e-10);

    return ok;
}