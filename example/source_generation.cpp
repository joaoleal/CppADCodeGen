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
#include <iosfwd>
#include <vector>
#include <cppadcg/cg.hpp>

using namespace CppAD;

int main(void) {
    // use a special object for source code generation
    typedef CG<double> CGD;
    typedef AD<CGD> ADCG;

    // independent variable vector
    CppAD::vector<ADCG> U(2);
    U[0] = 2.;
    U[1] = 3.;
    Independent(U);

    // dependent variable vector 
    CppAD::vector<ADCG> Z(1);

    // the model
    ADCG a = U[0] / 1. + U[1] * U[1];
    Z[0] = a / 2;

    ADFun<CGD> fun(U, Z); // the model tape

    // independent variable vector values
    CppAD::vector<double> u(2);
    u[0] = 2.;
    u[1] = 3.;

    /**
     * start the special steps for source code generation
     * for a Jacobian
     */
    CodeHandler<double> handler;

    CppAD::vector<CGD> indVars(2);
    handler.makeVariables(indVars);

    CppAD::vector<CGD> jac = fun.SparseJacobian(indVars);

    CLanguage<double> langC("double");
    CLangDefaultVariableNameGenerator<double> nameGen;

    std::ostringstream code;
    handler.generateCode(code, langC, jac, nameGen);
    std::cout << code.str();
}