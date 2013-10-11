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
#include "pow.hpp"

using namespace CppAD;

TEST_F(CppADCGOperationTest, PowTestOne) {
    // domain space vector
    size_t n = 2;
    double x = 0.5;
    double y = 2.;
    std::vector<double> u(n);
    u[0] = x;
    u[1] = y;

    test0nJac("PowTestOne", &PowTestOneFunc<double >, &PowTestOneFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, PowTestTwo) {
    // independent variable vector, indices, values, and declaration
    std::vector<double> u(2);
    size_t s = 0;
    size_t t = 1;
    u[s] = 2.;
    u[t] = 3.;

    test0nJac("PowTestTwo", &PowTestTwoFunc<double >, &PowTestTwoFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, PowTestThree) {
    // domain space vector
    size_t n = 1;
    std::vector<double> u(n);
    u[0] = 2.;

    test0nJac("PowTestThree", &PowTestThreeFunc<double >, &PowTestThreeFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, PowTestFour) {
    // domain space vector
    size_t n = 1;
    std::vector<double> u(n);
    u[0] = -2;

    test0nJac("PowTestFour", &PowTestFourFunc<double >, &PowTestFourFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, PowTestFive) {
    // domain space vector
    size_t n = 1;
    std::vector<double> u(n);
    u[0] = -1.;

    test0nJac("PowTestFive", &PowTestFiveFunc<double >, &PowTestFiveFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, PowTestSix) {
    // domain space vector
    size_t n = 1;
    std::vector<double> u(n);
    u[0] = 1.5;

    test0nJac("PowTestSix", &PowTestSixFunc<double >, &PowTestSixFunc<CG<double> >, u);
}
