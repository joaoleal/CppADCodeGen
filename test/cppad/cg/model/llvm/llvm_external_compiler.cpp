/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2019 Joao Leal
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

#include "LlvmModelTest.hpp"

using namespace CppAD;
using namespace CppAD::cg;

class LlvmModelExternalCompilerTest : public LlvmModelTest {
public:
    std::unique_ptr<LlvmModelLibrary<Base> > compileLib(LlvmModelLibraryProcessor<double>& p) override {
        ClangCompiler<double> clang("/usr/bin/clang-" + p.getVersion());
        return p.create(clang);
    }
};


TEST_F(LlvmModelExternalCompilerTest, ForwardZero) {
    testForwardZeroResults(*model, *fun, x);
}

TEST_F(LlvmModelExternalCompilerTest, DenseJacobian) {
    testDenseJacResults(*model, *fun, x);
}

TEST_F(LlvmModelExternalCompilerTest, DenseHessian) {
    testDenseHessianResults(*model, *fun, x);
}

TEST_F(LlvmModelExternalCompilerTest, Jacobian) {
    testJacobianResults(*llvmModelLib, *model, *fun, x, false);
}

TEST_F(LlvmModelExternalCompilerTest, Hessian) {
    testHessianResults(*llvmModelLib, *model, *fun, x, false);
}
