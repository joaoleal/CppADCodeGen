/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppadcg/cg.hpp>

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

bool Dynamic1() {
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
    CLangCompileModelHelper<double> compHelp(&fun, "dynamic");

    compHelp.setCreateForwardZero(true);
    compHelp.setCreateJacobian(true);
    compHelp.setCreateHessian(true);
    compHelp.setCreateSparseJacobian(true);
    compHelp.setCreateSparseHessian(true);
    compHelp.setMaxAssignmentsPerFunc(1);

    GccCompiler<double> compiler;

    CLangCompileDynamicHelper<double> compDynHelp(&compHelp);
    DynamicLib<double>* dynamicLib = compDynHelp.createDynamicLibrary(compiler);

    /**
     * test the library
     */
    bool ok = true;

    DynamicLibModel<double>* model = dynamicLib->model("dynamic");
    if (model == NULL)
        return false;

    // dimensions
    ok &= model->Domain() == fun.Domain();
    ok &= model->Range() == fun.Range();

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

    std::vector<double> depCGen = model->ForwardZero(x);
    ok &= compareValues("ForwardZero", depCGen, dep);

    // Jacobian
    std::vector<CGD> jac = fun.Jacobian(x2);
    depCGen = model->Jacobian(x);
    ok &= compareValues("Jacobian", depCGen, jac);

    // Hessian
    std::vector<CGD> w2(Z.size(), 1.0);
    std::vector<double> w(Z.size(), 1.0);

    std::vector<CGD> hess = fun.Hessian(x2, w2);
    depCGen = model->Hessian(x, w);
    ok &= compareValues("Hessian", depCGen, hess);

    // sparse Jacobian
    std::vector<double> jacCGen;
    std::vector<size_t> row, col;
    model->SparseJacobian(x, jacCGen, row, col);
    std::vector<double> jacCGenDense(jac.size());
    for (size_t i = 0; i < jacCGen.size(); i++) {
        jacCGenDense[row[i] * x.size() + col[i]] = jacCGen[i];
    }

    ok &= compareValues("sparse Jacobian", jacCGenDense, jac);

    // sparse Hessian
    std::vector<double> hessCGen;
    model->SparseHessian(x, w, hessCGen, row, col);
    std::vector<double> hessCGenDense(hess.size());
    for (size_t i = 0; i < hessCGen.size(); i++) {
        hessCGenDense[row[i] * x.size() + col[i]] = hessCGen[i];
    }

    ok &= compareValues("sparse Hessian", hessCGenDense, hess);

    delete dynamicLib;

    return ok;
}

bool Dynamic2() {
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
    CLangCompileModelHelper<double> compHelp(&fun, "dynamic2");

    compHelp.setCreateSparseJacobian(true);
    std::vector<size_t> row(3), col(3); // all elements except 1
    row[0] = 0;
    col[0] = 0;
    row[1] = 1;
    col[1] = 0;
    row[2] = 1;
    col[2] = 2;
    compHelp.setCustomSparseJacobianElements(row, col);

    compHelp.setCreateSparseHessian(true);
    row.resize(2);
    col.resize(2); // all elements except 1
    row[0] = 0;
    col[0] = 0;
    row[1] = 2;
    col[1] = 1;
    compHelp.setCustomSparseHessianElements(row, col);

    GccCompiler<double> compiler;
    compiler.setSourcesFolder("cppadcg_sources_2");

    CLangCompileDynamicHelper<double> compDynHelp(&compHelp);
    compDynHelp.setLibraryName("cppad_cg_model_2.so");
    DynamicLib<double>* dynamicLib = compDynHelp.createDynamicLibrary(compiler);

    /**
     * test the library
     */
    DynamicLibModel<double>* model = dynamicLib->model("dynamic2");
    if (model == NULL)
        return false;

    bool ok = true;

    // dimensions
    ok &= model->Domain() == fun.Domain();
    ok &= model->Range() == fun.Range();

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

    // sparse Jacobian
    std::vector<double> jacCGen;
    model->SparseJacobian(x, jacCGen, row, col);
    std::vector<CG<double> > jacSparse(row.size());

    std::vector<CGD> jac = fun.Jacobian(x2);
    for (size_t i = 0; i < row.size(); i++) {
        jacSparse[i] = jac[row[i] * x.size() + col[i]];
    }

    ok &= compareValues("sparse Jacobian", jacCGen, jacSparse);

    // sparse Hessian
    std::vector<double> w(Z.size(), 1.0);
    std::vector<double> hessCGen;
    model->SparseHessian(x, w, hessCGen, row, col);
    std::vector<CG<double> > hessSparse(row.size());

    std::vector<CGD> w2(Z.size(), 1.0);
    std::vector<CGD> hess = fun.Hessian(x2, w2);
    for (size_t i = 0; i < row.size(); i++) {
        hessSparse[i] = hess[row[i] * x.size() + col[i]];
    }

    ok &= compareValues("sparse Hessian", hessCGen, hessSparse);

    delete dynamicLib;

    return ok;
}

bool Dynamic() {
    bool ok = true;
    ok &= Dynamic1();
    ok &= Dynamic2();

    return ok;
}