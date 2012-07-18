/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
#include <iosfwd>
#include <vector>
#include <cppad_cgoo/cg.hpp>

using namespace CppAD;

int main(void) {
    // use a special object for source code generation
    typedef CG<double> CGD;
    typedef AD<CGD> ADCG;

    // independent variable vector
    std::vector<ADCG> U(2);
    U[0] = 2.;
    U[1] = 3.;
    Independent(U);

    // dependent variable vector 
    std::vector<ADCG> Z(1);

    // the model
    ADCG a = U[0] / 1. + U[1] * U[1];
    Z[0] = a / 2;

    ADFun<CGD> fun(U, Z);

    // independent variable vector, indices, values, and declaration
    std::vector<double> u(2);
    u[0] = 2.;
    u[1] = 3.;

    /**
     * start the special steps for source code generation
     * for a jacobian
     */
    CodeHandler<double> handler;

    std::vector<CGD> indVars(2);
    handler.makeVariables(indVars);

    std::vector<CGD> jac = fun.SparseJacobian(indVars);

    CLanguage<double> langC("double");
    CLangDefaultVariableNameGenerator<double> nameGen;

    std::ostringstream code;
    handler.generateCode(code, langC, jac, nameGen);
    std::cout << code.str();
}