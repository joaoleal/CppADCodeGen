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
#include <cppadcg/cg.hpp>

using namespace CppAD;

int main(void) {
    // use a special object for source code generation
    typedef CG<double> CGD;
    typedef AD<CGD> ADCG;

    // independent variable vector
    std::vector<ADCG> U(2);
    Independent(U);

    // dependent variable vector 
    std::vector<ADCG> Z(1);

    // the model
    ADCG a = U[0] / 1. + U[1] * U[1];
    Z[0] = a / 2;

    ADFun<CGD> fun(U, Z);

    /**
     * Create the dynamic library
     * (generates and compiles source code)
     */
    CLangCompileModelHelper<double> compModelH(&fun, "model");
    compModelH.setCreateJacobian(true);

    CLangCompileDynamicHelper<double> compDynH(&compModelH);

    GccCompiler<double> compiler;
    DynamicLib<double>* dynamicLib = compDynH.createDynamicLibrary(compiler);

    /**
     * Use the dynamic library
     */
    DynamicLibModel<double>* model = dynamicLib->model("model");
    std::vector<double> x(U.size());
    x[0] = 2.5;
    x[1] = 3.5;
    std::vector<double> jac = model->Jacobian(x);

    // print out the result
    std::cout << jac[0] << " " << jac[1] << std::endl;

    delete model;
    delete dynamicLib;
}
