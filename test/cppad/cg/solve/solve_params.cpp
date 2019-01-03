/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2019 Joao Leal
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */
#include "CppADCGSolveTest.hpp"

using namespace CppAD;
using namespace CppAD::cg;

TEST_F(CppADCGSolveTest, SolveParams) {
    // independent variable vector
    std::vector<ADCGD> x(2);
    x[0] = -1.5;
    x[1] = -0.5;

    std::vector<ADCGD> p(2);
    p[0] = 1.0;
    p[1] = 3.0;

    // use a special object for source code generation
    // declare independent variables, dynamic parameters, starting recording
    size_t abort_op_index = 0;
    bool record_compare = true;
    CppAD::Independent(x, abort_op_index, record_compare, p);

    // dependent variable vector
    std::vector<ADCGD> y(3);

    // model
    y[0] = x[0] + x[1];
    y[1] = y[0] + p[0];
    y[2] = p[1] + y[1] - 2;

    // create f: x -> y
    ADFun<CGD> fun(x, y);

    test_solve(fun, 2, 1, x, p);
}
