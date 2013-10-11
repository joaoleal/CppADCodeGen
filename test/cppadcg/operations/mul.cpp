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
#include "mul.hpp"

using namespace CppAD;

TEST_F(CppADCGOperationTest, MulTestOne) {
    // independent variable vector, indices, values, and declaration
    std::vector<double> u(2);
    size_t s = 0;
    size_t t = 1;
    u[s] = 3.;
    u[t] = 2.;

    test0nJac("MulTestOne", &MulTestOneFunc<double >, &MulTestOneFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, MulTestTwo) {
    // independent variable vector
    std::vector<double> u(1);
    u[0] = .5;

    test0nJac("MulTestTwo", &MulTestTwoFunc<double >, &MulTestTwoFunc<CG<double> >, u);
}
