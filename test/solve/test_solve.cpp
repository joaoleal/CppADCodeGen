/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppadcg/cg.hpp>
#include <assert.h>
#include <limits>

#include "test_solve.hpp"

using namespace CppAD;
using namespace std;

typedef CG<double> CGD;

extern bool test_verbose;
extern bool test_printvalues;

bool test_solve(ADFun<CGD>& fun,
                size_t expressionIndex,
                size_t indIndex,
                const std::vector<AD<CGD> >& testValues) {

    std::vector<double> testValuesD(testValues.size());
    for (size_t i = 0; i < testValues.size(); i++) {
        testValuesD[i] = CppAD::Value(CppAD::Var2Par(testValues[i])).getParameterValue();
    }

    return test_solve(fun, expressionIndex, indIndex, testValuesD);
}

bool test_solve(ADFun<CGD>& fun,
                size_t expressionIndex,
                size_t indIndex,
                const std::vector<double>& testValues) {

    using std::vector;

    size_t m = fun.Range();
    size_t n = fun.Domain();

    assert(expressionIndex < m);
    assert(indIndex < n);
    assert(testValues.size() == n);

    // evaluate the dependent values
    const vector<double> depValues = calculateDependentForward0(fun, testValues);

    assert(NearEqual(depValues[expressionIndex], 0.0, 10e-4, 10e-4));

    // generate the operation graph
    CodeHandler<double> handler;
    vector<CGD> indVars(n);
    handler.makeVariables(indVars);

    vector<CGD> dep = fun.Forward(0, indVars);

    /**
     * solve
     */
    CGD solution = handler.solveFor(dep[expressionIndex].getSourceCodeFragment(), indVars[indIndex].getSourceCodeFragment());

    if (test_verbose) {
        printModel(handler, solution);
    }

    /**
     * create a new tape
     */
    vector<CGD> newDep(1);
    newDep[0] = solution;

    // new independent vector (without one variable)
    vector<AD<double> > newIndep(n - 1);
    for (size_t i = 0; i < indIndex; i++) {
        newIndep[i] = testValues[i];
    }
    for (size_t i = indIndex + 1; i < n; i++) {
        newIndep[i - 1] = testValues[i];
    }

    // create a longer independent vector to use in the evaluator
    vector<AD<double> > newIndepLong(n);
    for (size_t i = 0; i < indIndex; i++) {
        newIndepLong[i] = newIndep[i];
    }
    newIndepLong[indIndex] = std::numeric_limits<double>::quiet_NaN();
    for (size_t i = indIndex + 1; i < n; i++) {
        newIndepLong[i] = newIndep[i - 1];
    }

    // determine the result
    Evaluator<double, double> evaluator(handler, newDep);
    vector<AD<double> > result = evaluator.evaluate(newIndepLong);

    double resultVal = CppAD::Value(CppAD::Var2Par(result[0]));

    return NearEqual(resultVal, testValues[indIndex], 10e-4, 10e-4);
}

std::vector<double> calculateDependentForward0(ADFun<CGD>& fun, const std::vector<double>& testValues) {
    using std::vector;

    size_t m = fun.Range();
    size_t n = fun.Domain();

    vector<CGD> indVars(n);
    for (size_t i = 0; i < n; ++i) {
        indVars[i] = testValues[i];
    }

    vector<CGD> dep = fun.Forward(0, indVars);

    vector<double> depVals(m);
    for (size_t i = 0; i < m; ++i) {
        depVals[i] = dep[i].getParameterValue();
    }

    return depVals;
}

void printModel(CodeHandler<double>&handler, CGD& dep) {
    CLanguage<double> langC("double");
    CLangDefaultVariableNameGenerator<double> nameGen;

    std::vector<CGD> depv(1);
    depv[0] = dep;

    std::ostringstream code;
    handler.generateCode(code, langC, depv, nameGen);
    std::cout << code.str();
}