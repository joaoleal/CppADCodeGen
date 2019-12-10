#ifndef CPPAD_CG_TEST_CPPADCGMODELTEST_INCLUDED
#define CPPAD_CG_TEST_CPPADCGMODELTEST_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2019 Joao Leal
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
#include "CppADCGTest.hpp"

namespace CppAD {
namespace cg {

class CppADCGModelTest : public CppADCGTest {
public:
    using CGD = CG<double>;
    using ADCG = AD<CGD>;
public:

    explicit CppADCGModelTest(bool verbose = false,
                              bool printValues = false) :
            CppADCGTest(verbose, printValues) {
    }

    static std::vector<CGD> makeVector(const std::vector<Base>& x) {
        std::vector<CGD> x2(x.size());
        for (size_t i = 0; i < x.size(); i++) {
            x2[i] = x[i];
        }
        return x2;
    }

    /**
     * Compares the results from forward zero.
     *
     * @param model
     * @param fun
     * @param x independent vector values
     * @param epsilonR relative error
     * @param epsilonA absolute error
     */
    void testForwardZeroResults(GenericModel<Base>& model,
                                ADFun<CGD>& fun,
                                const std::vector<Base>& x,
                                double epsilonR = 1e-14,
                                double epsilonA = 1e-14) {
        ASSERT_EQ(model.Domain(), fun.Domain());
        ASSERT_EQ(model.Range(), fun.Range());

        std::vector<CGD> x2 = makeVector(x);

        std::vector<CGD> dep = fun.Forward(0, x2);
        std::vector<Base> depCGen = model.ForwardZero(x);

        ASSERT_TRUE(compareValues(depCGen, dep, epsilonR, epsilonA));
    }

    /**
     * Compares the results from a dense Jacobian.
     *
     * @param model
     * @param fun
     * @param x independent vector values
     * @param epsilonR relative error
     * @param epsilonA absolute error
     */
    void testDenseJacResults(GenericModel<Base>& model,
                             ADFun<CGD>& fun,
                             const std::vector<Base>& x,
                             double epsilonR = 1e-14,
                             double epsilonA = 1e-14) {
        ASSERT_EQ(model.Domain(), fun.Domain());
        ASSERT_EQ(model.Range(), fun.Range());

        std::vector<CGD> x2 = makeVector(x);

        auto jac = fun.Jacobian(x2);
        if (verbose_) {
            std::cout << "ADFun Jacobian" << std::endl;
            print(jac);
        }

        auto depCGen = model.Jacobian(x);

        ASSERT_TRUE(compareValues(depCGen, jac, epsilonR, epsilonA));
    }

    /**
     * Compares the results from a Dense Hessian.
     *
     * @param model
     * @param fun
     * @param x independent vector values
     * @param epsilonR relative error
     * @param epsilonA absolute error
     */
    void testDenseHessianResults(GenericModel<Base>& model,
                                 ADFun<CGD>& fun,
                                 const std::vector<Base>& x,
                                 double epsilonR = 1e-14,
                                 double epsilonA = 1e-14) {
        // dimensions
        ASSERT_EQ(model.Domain(), fun.Domain());
        ASSERT_EQ(model.Range(), fun.Range());

        std::vector<CGD> x2 = makeVector(x);
        std::vector<CGD> w2(fun.Range(), 1.0);
        std::vector<Base> w(fun.Range(), 1.0);

        std::vector<CGD> hess = fun.Hessian(x2, w2);

        if (verbose_) {
            std::cout << "ADFun Hessian" << std::endl;
            print(hess);
        }

        auto depCGen = model.Hessian(x, w);

        ASSERT_TRUE(compareValues(depCGen, hess, epsilonR, epsilonA));
    }

private:
    void testJacobianResults(GenericModel<Base>& model,
                             const std::vector<CGD>& jac,
                             const std::vector<Base>& x,
                             bool customSparsity,
                             double epsilonR = 1e-14,
                             double epsilonA = 1e-14) {
        std::vector<Base> jacCGen;
        std::vector<size_t> row, col;

        model.SparseJacobian(x, jacCGen, row, col);

        std::vector<Base> jacCGenDense(jac.size());
        
        for (size_t i = 0; i < jacCGen.size(); i++) {
            size_t p = row[i] * x.size() + col[i];
            jacCGenDense[p] = jacCGen[i];
        }

        std::vector<CGD> jacAdFunPartial(jac.size());
        if (customSparsity) {
            for (size_t i = 0; i < jacCGen.size(); i++) {
                size_t p = row[i] * x.size() + col[i];
                jacAdFunPartial[p] = jac[p];
            }
        } else {
            jacAdFunPartial = jac;
        }

        if (verbose_) {
            std::cout << "sparse Jacobian (in a dense format)" << std::endl;
            print(jacCGenDense);
        }

        ASSERT_TRUE(this->compareValues(jacCGenDense, jacAdFunPartial, epsilonR, epsilonA));
    }

    void testHessianResults(GenericModel<Base>& model,
                            const std::vector<CGD>& hess,
                            const std::vector<Base>& x,
                            bool customSparsity,
                            double epsilonR = 1e-14,
                            double epsilonA = 1e-14) {
        std::vector<CGD> w2(model.Range(), 1.0);
        std::vector<Base> w(model.Range(), 1.0);
        std::vector<size_t> row, col;
        std::vector<Base> hessCGen;

        model.SparseHessian(x, w, hessCGen, row, col);

        std::vector<Base> hessCGenDense(hess.size());
        for (size_t i = 0; i < hessCGen.size(); i++) {
            size_t p = row[i] * x.size() + col[i];
            hessCGenDense[p] = hessCGen[i];
        }

        std::vector<CGD> hessAdFunPartial(hess.size());
        if (customSparsity) {
            for (size_t i = 0; i < hessCGen.size(); i++) {
                size_t p = row[i] * x.size() + col[i];
                hessAdFunPartial[p] = hess[p];
            }
        } else {
            hessAdFunPartial = hess;
        }

        if (verbose_) {
            std::cout << "sparse Hessian (in a dense format)" << std::endl;
            print(hessCGenDense);
        }

        ASSERT_TRUE(this->compareValues(hessCGenDense, hessAdFunPartial, epsilonR, epsilonA));
    }

public:

    /**
     * Compares the results from a sparse Jacobian.
     *
     * @param model
     * @param fun
     * @param x independent vector values
     * @param epsilonR relative error
     * @param epsilonA absolute error
     */
    void testJacobianResults(ModelLibrary<Base>& lib,
                             GenericModel<Base>& model,
                             ADFun<CGD>& fun,
                             const std::vector<Base>& x,
                             bool customSparsity,
                             double epsilonR = 1e-14,
                             double epsilonA = 1e-14) {
        ASSERT_EQ(model.Domain(), fun.Domain());
        ASSERT_EQ(model.Range(), fun.Range());

        std::vector<CGD> x2 = makeVector(x);
        std::vector<CGD> jac = fun.Jacobian(x2);

        if (verbose_) {
            std::cout << "ADFun Jacobian" << std::endl;
            print(jac);
        }

        testJacobianResults(model, jac, x, customSparsity, epsilonR, epsilonA);

        if (lib.getThreadNumber() > 1) {
            // sparse Jacobian again (make sure the second run is also OK)
            testJacobianResults(model, jac, x, customSparsity, epsilonR, epsilonA);
        }

    }

    /**
     * Compares the results from a sparse Hessian.
     *
     * @param model
     * @param fun
     * @param x independent vector values
     * @param epsilonR relative error
     * @param epsilonA absolute error
     */
    void testHessianResults(ModelLibrary<Base>& lib,
                            GenericModel<Base>& model,
                            ADFun<CGD>& fun,
                            const std::vector<Base>& x,
                            bool customSparsity,
                            double epsilonR = 1e-14,
                            double epsilonA = 1e-14) {
        ASSERT_EQ(model.Domain(), fun.Domain());
        ASSERT_EQ(model.Range(), fun.Range());

        std::vector<CGD> x2 = makeVector(x);
        std::vector<CGD> w2(fun.Range(), 1.0);
        std::vector<Base> w(fun.Range(), 1.0);

        std::vector<CGD> hess = fun.Hessian(x2, w2);

        if (verbose_) {
            std::cout << "ADFun Hessian" << std::endl;
            print(hess);
        }

        testHessianResults(model, hess, x, customSparsity, epsilonR, epsilonA);

        if (lib.getThreadNumber() > 1) {
            // sparse Hessian again (make sure the second run is also OK)
            testHessianResults(model, hess, x, customSparsity, epsilonR, epsilonA);
        }
    }

    /**
     * Compares the results from Hessian, Jacobian, sparse Hessian and 
     * sparse Jacobian.
     * 
     * @param model
     * @param fun
     * @param x independent vector values
     * @param epsilonR relative error
     * @param epsilonA absolute error
     */
    void testModelResults(ModelLibrary<Base>& lib,
                          GenericModel<Base>& model,
                          ADFun<CGD>& fun,
                          const std::vector<Base>& x,
                          bool customSparsity,
                          double epsilonR = 1e-14,
                          double epsilonA = 1e-14,
                          bool denseJacobian = true,
                          bool denseHessian = true) {
        // dimensions
        ASSERT_EQ(model.Domain(), fun.Domain());
        ASSERT_EQ(model.Range(), fun.Range());

        testForwardZeroResults(model, fun, x, epsilonR, epsilonA);

        // Jacobian
        if (denseJacobian) {
            testDenseJacResults(model, fun, x, epsilonR, epsilonA);
        }

        if (denseHessian) {
            testDenseHessianResults(model, fun, x, epsilonR, epsilonA);
        }

        // sparse Jacobian
        testJacobianResults(lib, model, fun, x, customSparsity, epsilonR, epsilonA);

        // sparse Hessian
        testHessianResults(lib, model, fun, x, customSparsity, epsilonR, epsilonA);
    }

    inline ::testing::AssertionResult compareValues(const std::vector<double>& depCGen,
                                                    const std::vector<CppAD::cg::CG<double> >& dep,
                                                    double epsilonR = 1e-14, double epsilonA = 1e-14) {

        std::vector<double> depd(dep.size());

        for (size_t i = 0; i < depd.size(); i++) {
            depd[i] = dep[i].getValue();
        }

        return CppADCGTest::compareValues(depCGen, depd, epsilonR, epsilonA);
    }

};

} // END cg namespace
} // END CppAD namespace

#endif