/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad_cgoo/cg.hpp>

#include "gcc_load_dynamic.hpp"

using namespace CppAD;
using namespace std;

bool compareValues(const std::string& testType,
                   const std::vector<double>& depCGen,
                   const std::vector<CG<double> >& dep,
                   double epsilonR = 1e-14, double epsilonA = 1e-14) {

    std::vector<double> depd(dep.size());

    for (size_t i = 0; i < depd.size(); i++) {
        depd[i] = dep[i].getParameterValue();
    }

    return compareValues(testType, depCGen, depd, epsilonR, epsilonA);
}

bool Dynamic() {
    // use a special object for source code generation
    typedef CG<double> CGD;
    typedef AD<CGD> ADCG;

    // independent variables
    std::vector<ADCG> u(3);
    u[0] = 1;
    u[1] = 1;
    u[2] = 1;

    CppAD::Independent(u);

    // dependent variable vector 
    std::vector<ADCG> Z(2);

    /**
     * create the CppAD tape as usual
     */
    Z[0] = cos(u[0]);
    Z[1] = u[1] * u[2] + sin(u[0]);

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD> fun(u, Z);

    /**
     * Create the dynamic library
     * (generate and compile source code)
     */
    CLangCompileHelper<double> compHelp(&fun);

    compHelp.setCreateForwardZero(true);
    compHelp.setCreateJacobian(true);
    compHelp.setCreateHessian(true);
    compHelp.setCreateSparseJacobian(true);
    compHelp.setCreateSparseHessian(true);

    GccCompiler<double> compiler;

    DynamicLib<double>* dynamicLib = compHelp.createDynamicLibrary(compiler);

    /**
     * test the library
     */
    bool ok = true;

    // dimensions
    ok &= dynamicLib->Domain() == fun.Domain();
    ok &= dynamicLib->Range() == fun.Range();

    /**
     */
    std::vector<double> x(u.size());
    x[0] = 1;
    x[1] = 2;
    x[2] = 1;

    std::vector<CGD> x2(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        x2[i] = x[i];
    }

    // forward zero
    std::vector<CGD> dep = fun.Forward(0, x2);

    std::vector<double> depCGen = dynamicLib->ForwardZero(x);
    ok &= compareValues("ForwardZero", depCGen, dep);

    // Jacobian
    std::vector<CGD> jac = fun.Jacobian(x2);
    depCGen = dynamicLib->Jacobian(x);
    ok &= compareValues("Jacobian", depCGen, jac);

    // Hessian
    std::vector<CGD> w2(Z.size(), 1.0);
    std::vector<double> w(Z.size(), 1.0);

    std::vector<CGD> hess = fun.Hessian(x2, w2);
    depCGen = dynamicLib->Hessian(x, w);
    ok &= compareValues("Hessian", depCGen, hess);

    // sparse Jacobian
    std::vector<double> jacCGen;
    std::vector<size_t> row, col;
    dynamicLib->SparseJacobian(x, jacCGen, row, col);
    std::vector<double> jacCGenDense(jac.size());
    for (size_t i = 0; i < jacCGen.size(); i++) {
        jacCGenDense[row[i] * x.size() + col[i]] = jacCGen[i];
    }

    ok &= compareValues("sparse Jacobian", jacCGenDense, jac);

    // sparse Hessian
    std::vector<double> hessCGen;
    dynamicLib->SparseHessian(x, w, hessCGen, row, col);
    std::vector<double> hessCGenDense(hess.size());
    for (size_t i = 0; i < hessCGen.size(); i++) {
        hessCGenDense[row[i] * x.size() + col[i]] = hessCGen[i];
    }

    ok &= compareValues("sparse Hessian", hessCGenDense, hess);

    delete dynamicLib;

    return ok;
}