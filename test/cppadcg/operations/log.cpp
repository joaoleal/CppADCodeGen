/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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
#include "CppADCGOperationTest.hpp"
#include "log.hpp"

using namespace CppAD;

TEST_F(CppADCGOperationTest, LogTestOne) {
    // independent variable vector, indices, values, and declaration
    std::vector<double> u(1);
    u[0] = 2.;

    test0nJac("LogTestOne", &LogTestOneFunc<double >, &LogTestOneFunc<CG<double> >, u, 1e-10, 1e-10);
}

TEST_F(CppADCGOperationTest, LogTestTwo) {
    // independent variable vector
    std::vector<double> u(1);
    u[0] = 1.;

    test0nJac("LogTestTwo", &LogTestTwoFunc<double >, &LogTestTwoFunc<CG<double> >, u, 1e-10, 1e-10);
}