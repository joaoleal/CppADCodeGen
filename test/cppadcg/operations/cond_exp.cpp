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
#include "cond_exp.hpp"

using namespace CppAD;

TEST_F(CppADCGOperationTest, CondExp_pvvv) {
    // independent variable vector
    std::vector<double> u0(3);
    u0[0] = 0.;
    u0[1] = 1.;
    u0[2] = 2.;
    std::vector<double> u1(3);
    u1[0] = 1.5;
    u1[1] = 0.;
    u1[2] = 2.;

    std::vector<std::vector<double> > u(2);
    u[0] = u0;
    u[1] = u1;

    test0nJac("CondExp_pvvv", &CondExp_pvvvFunc<double >, &CondExp_pvvvFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, CondExp_vpvv) {
    // independent variable vector
    std::vector<double> u(3);
    u[0] = 0.;
    u[1] = 1.;
    u[2] = 2.;

    test0nJac("CondExp_vpvv", &CondExp_vpvvFunc<double >, &CondExp_vpvvFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, CondExp_vvpv) {
    // independent variable vector
    std::vector<double> u(3);
    u[0] = 0.;
    u[1] = 1.;
    u[2] = 2.;

    test0nJac("CondExp_vvpv", &CondExp_vvpvFunc<double >, &CondExp_vvpvFunc<CG<double> >, u);
}

TEST_F(CppADCGOperationTest, CondExp_vvvp) {
    // independent variable vector
    std::vector<double> u(3);
    u[0] = 0.;
    u[1] = 1.;
    u[2] = 2.;

    test0nJac("CondExp_vvvp", &CondExp_vvvpFunc<double >, &CondExp_vvvpFunc<CG<double> >, u);
}
