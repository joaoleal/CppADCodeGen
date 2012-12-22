#ifndef TEST_SOLVE_INCLUDED
#define	TEST_SOLVE_INCLUDED
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

bool test_solve(CppAD::ADFun<CppAD::CG<double> >& fun,
                size_t expressionIndex,
                size_t indIndex,
                const std::vector<CppAD::AD<CppAD::CG<double> > >& testValues);

bool test_solve(CppAD::ADFun<CppAD::CG<double> >& fun,
                size_t expressionIndex,
                size_t indIndex,
                const std::vector<double>& testValues);

std::vector<double> calculateDependentForward0(CppAD::ADFun<CppAD::CG<double> >& fun,
                                               const std::vector<double>& testValues);

void printModel(CppAD::CodeHandler<double>& handler, CppAD::CG<double>& dep);

#endif
