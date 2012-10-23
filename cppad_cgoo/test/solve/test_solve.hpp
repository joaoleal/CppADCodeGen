#ifndef TEST_SOLVE_INCLUDED
#define	TEST_SOLVE_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */


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
