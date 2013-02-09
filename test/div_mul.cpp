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
#include "CppADCGOperationTest.hpp"
#include "div_mul.hpp"

using namespace CppAD;

TEST_F(CppADCGOperationTest, DivMulTestOne) {
    // independent variable vector, indices, values, and declaration
    std::vector<double> u(3);
    u[0] = 2.;
    u[1] = 3.;
    u[2] = 4.;

    test0nJac("DivMulTestOne", &DivMulTestOneFunc<double >, &DivMulTestOneFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, DivMulTestTwo) {
    // independent variable vector
    std::vector<double> u(4);
    u[0] = 2.;
    u[1] = 3.;
    u[2] = 4.;
    u[3] = 5.;

    test0nJac("DivMulTestTwo", &DivMulTestTwoFunc<double >, &DivMulTestTwoFunc<CG<double> >, u);
}
