/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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
#include <cppad/cg/cppadcg.hpp>
#include <cppad/cg/model/llvm/llvm.hpp>

#include "CppADCGModelTest.hpp"

using namespace CppAD;
using namespace CppAD::cg;

class LlvmModelTest : public CppAD::cg::CppADCGModelTest {
public:

    // Per-test-case tear-down.
    // Called after the last test in this test case.
    static void TearDownTestCase() {
        // llvm_shutdown must be called once for all tests!
        llvm::llvm_shutdown();
    }
};

template<class T>
CppAD::ADFun<T>* modelFunc(std::vector<CppAD::AD<T> >& x) {
    using namespace CppAD;
    using namespace std;

    assert(x.size() == 3);

    CppAD::Independent(x);

    std::vector<CppAD::AD<T> > y(2);
    y[0] = CppAD::abs(x[0]) * x[1];
    y[1] = CppAD::cos(x[2]);

    // f(v) = |w|
    return new CppAD::ADFun<T>(x, y);
}

TEST_F(LlvmModelTest, llvm) {
    std::vector<double> x(3);
    x[0] = -1;
    x[1] = 2;
    x[2] = 3;

    std::vector<AD<CG<double> > > u(3);
    //u[0] = x[0];

    std::unique_ptr<CppAD::ADFun<CG<Base> > > fun(modelFunc<CG<Base> >(u));

    /**
     * Create the dynamic library
     * (generate and compile source code)
     */
    ModelCSourceGen<double> compHelp(*fun.get(), "mySmallModel");
    compHelp.setCreateForwardZero(true);
    compHelp.setCreateJacobian(true);
    compHelp.setCreateHessian(true);
    compHelp.setCreateSparseJacobian(true);
    compHelp.setCreateSparseHessian(true);
    compHelp.setCreateForwardOne(true);
    compHelp.setMultiThreading(false);

    ModelLibraryCSourceGen<double> compDynHelp(compHelp);
    compDynHelp.setVerbose(this->verbose_);
    compDynHelp.setMultiThreading(MultiThreadingType::NONE);

    LlvmModelLibraryProcessor<double> p(compDynHelp);

    std::unique_ptr<LlvmModelLibrary<Base> > llvmModelLib(p.create());
    std::unique_ptr<GenericModel<Base> > model(llvmModelLib->model("mySmallModel"));
    ASSERT_TRUE(model.get() != nullptr);

    this->testModelResults(*llvmModelLib, *model, *fun.get(), x);

    model.reset(nullptr); // must be freed before llvm_shutdown()
    llvmModelLib.reset(nullptr); // must be freed before llvm_shutdown()
}

TEST_F(LlvmModelTest, llvm_externalCompiler) {
    std::vector<double> x(3);
    x[0] = -1;
    x[1] = 2;
    x[2] = 3;

    std::vector<AD<CG<double> > > u(3);
    //u[0] = x[0];

    std::unique_ptr<CppAD::ADFun<CG<Base> > > fun(modelFunc<CG<Base> >(u));

    /**
     * Create the dynamic library
     * (generate and compile source code)
     */
    ModelCSourceGen<double> compHelp(*fun.get(), "mySmallModel");
    compHelp.setCreateForwardZero(true);
    compHelp.setCreateJacobian(true);
    compHelp.setCreateHessian(true);
    compHelp.setCreateSparseJacobian(true);
    compHelp.setCreateSparseHessian(true);
    compHelp.setCreateForwardOne(true);
    compHelp.setMultiThreading(false);

    ModelLibraryCSourceGen<double> compDynHelp(compHelp);
    compDynHelp.setVerbose(this->verbose_);
    compDynHelp.setMultiThreading(MultiThreadingType::NONE);

    LlvmModelLibraryProcessor<double> p(compDynHelp);

    ClangCompiler<double> clang;

    std::unique_ptr<LlvmModelLibrary<Base> > llvmModelLib(p.create(clang));
    std::unique_ptr<GenericModel<Base> > model(llvmModelLib->model("mySmallModel"));
    ASSERT_TRUE(model.get() != nullptr);

    this->testModelResults(*llvmModelLib, *model, *fun.get(), x);

    model.reset(nullptr); // must be freed before llvm_shutdown()
    llvmModelLib.reset(nullptr); // must be freed before llvm_shutdown()
}