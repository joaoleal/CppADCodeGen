#ifndef CPPAD_CG_TEST_CPPADCGDYNAMICATOMICTEST_INCLUDED
#define CPPAD_CG_TEST_CPPADCGDYNAMICATOMICTEST_INCLUDED
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

#include "CppADCGTest.hpp"
#include "gccCompilerFlags.hpp"

namespace CppAD {
namespace cg {

class CppADCGDynamicAtomicTest : public CppADCGTest {
public:
    typedef CppADCGTest::Base Base;
    typedef CppADCGTest::CGD CGD;
    typedef CppADCGTest::ADCGD ADCGD;
protected:
    const std::string _modelName;
    std::unique_ptr<ADFun<CGD>> _funInner;
    std::unique_ptr<ADFun<CGD>> _funOuter; // no atomics
    std::unique_ptr<DynamicLib<Base>> _dynamicLib;
    std::unique_ptr<DynamicLib<Base>> _dynamicLib2;
    std::unique_ptr<GenericModel<Base>> _modelLib;
public:

    inline CppADCGDynamicAtomicTest(const std::string& modelName,
                                    bool verbose = false,
                                    bool printValues = false) :
        CppADCGTest(verbose, printValues),
        _modelName(modelName),
        _funInner(nullptr),
        _funOuter(nullptr) {
        //this->verbose_ = true;
    }

    virtual std::vector<ADCGD> model(const std::vector<ADCGD>& x) = 0;

    virtual std::vector<ADCGD> modelOuter(const std::vector<ADCGD>& y) {
        std::vector<ADCGD> Z(y.size() - 1);

        for (size_t i = 0; i < y.size() - 1; i++) {
            Z[i] = 2 * y[i];
        }
        Z[Z.size() - 1] += y[y.size() - 1];

        return Z;
    }

    virtual void TearDown() {
        _dynamicLib.reset();
        _dynamicLib2.reset();
        _modelLib.reset();
        _funInner.reset();
        _funOuter.reset();
    }

    virtual ~CppADCGDynamicAtomicTest() {
    }

    void testADFunAtomicLib(const CppAD::vector<Base>& x) {
        CppAD::vector<Base> xNorm(x.size());
        for (size_t i = 0; i < xNorm.size(); i++)
            xNorm[i] = 1.0;
        CppAD::vector<Base> eqNorm;

        testADFunAtomicLib(x, xNorm, eqNorm);
    }

    /**
     * Tests one compiled model used as an atomic function by an ADFun
     * The outer model method (modelOuter) is NOT used here!
     * It uses z = g(x) = f(x).
     */
    void testADFunAtomicLib(const CppAD::vector<Base>& x,
                            const CppAD::vector<Base>& xNorm,
                            const CppAD::vector<Base>& eqNorm,
                            Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        ASSERT_EQ(x.size(), xNorm.size());

        using namespace CppAD;
        using namespace std;
        using CppAD::vector;

        prepareAtomicLib(x, xNorm, eqNorm);

        const size_t n = _funInner->Domain();
        const size_t m = _funInner->Range();

        vector< AD<double> > ax(n);
        for (size_t j = 0; j < n; j++)
            ax[j] = x[j];

        // declare independent variables and start tape recording
        CppAD::Independent(ax);

        vector<AD<double> > ay(m);

        // call user function and store CGAtomicLibModel(x) in au[0] 
        unique_ptr<GenericModel<Base> > modelLib(_dynamicLib->model(_modelName));
        CGAtomicGenericModel<double>& atomicfun = modelLib->asAtomic();

        atomicfun(ax, ay);

        // create f2: x -> y and stop tape recording
        ADFun<double> f2(ax, ay);

        /**
         * Test zero order
         */
        vector<CGD> xOrig(n);
        for (size_t j = 0; j < n; j++)
            xOrig[j] = x[j];

        vector<CGD> yOrig = _funInner->Forward(0, xOrig);
        vector<double> yInner = modelLib->ForwardZero(x);
        vector<double> yOuter = f2.Forward(0, x);

        ASSERT_TRUE(compareValues(yInner, yOrig, epsilonR, epsilonA));
        ASSERT_TRUE(compareValues(yOuter, yOrig, epsilonR, epsilonA));

        /**
         * Test first order forward mode
         */
        size_t k = 1;
        size_t k1 = k + 1;

        vector<double> x_p(n);
        for (size_t j = 0; j < n; j++)
            x_p[j] = 0;

        vector<CGD> x_pOrig(n);
        vector<double> tx(k1 * n);
        for (size_t j = 0; j < n; j++)
            tx[j * k1] = x[j]; // zero order
        for (size_t j = 0; j < n; j++)
            tx[j * k1 + 1] = 0; // first order

        for (size_t j = 0; j < n; j++) {
            x_p[j] = 1;
            x_pOrig[j] = 1;
            tx[j * k1 + 1] = 1;

            vector<CGD> y_pOrig = _funInner->Forward(1, x_pOrig);
            vector<double> y_pInner = modelLib->ForwardOne(tx);
            vector<double> y_pOuter = f2.Forward(1, x_p);

            x_p[j] = 0;
            x_pOrig[j] = 0;
            tx[j * k1 + 1] = 0;

            ASSERT_TRUE(compareValues(y_pInner, y_pOrig, epsilonR, epsilonA));
            ASSERT_TRUE(compareValues(y_pOuter, y_pOrig, epsilonR, epsilonA));
        }

        /**
         * Test first order reverse mode
         */
        k = 0;
        k1 = k + 1;

        vector<double> w(m);
        for (size_t i = 0; i < m; i++)
            w[i] = 0;
        vector<CGD> wOrig(m);
        tx.resize(k1 * n);
        for (size_t j = 0; j < n; j++)
            tx[j * k1] = x[j]; // zero order
        vector<double> ty(k1 * m);
        for (size_t i = 0; i < m; i++)
            ty[i * k1] = yInner[i]; // zero order

        for (size_t i = 0; i < m; i++) {
            w[i] = 1;
            wOrig[i] = 1;

            vector<CGD> dwOrig = _funInner->Reverse(1, wOrig);
            vector<double> dwInner = modelLib->ReverseOne(tx, ty, w);
            vector<double> dwOuter = f2.Reverse(1, w);

            w[i] = 0;
            wOrig[i] = 0;

            ASSERT_TRUE(compareValues(dwInner, dwOrig, epsilonR, epsilonA));
            ASSERT_TRUE(compareValues(dwOuter, dwOrig, epsilonR, epsilonA));
        }

        /**
         * Test second order reverse mode
         */
        k = 1;
        k1 = k + 1;
        tx.resize(k1 * n);
        ty.resize(k1 * m);
        vector<double> py(k1 * m);
        vector<CGD> pyOrig(k1 * m);
        //wOrig.resize(k1 * m);
        for (size_t j = 0; j < n; j++) {
            tx[j * k1] = x[j]; // zero order
            tx[j * k1 + 1] = 0; // first order
        }
        for (size_t i = 0; i < m; i++) {
            ty[i * k1] = yInner[i]; // zero order
            py[i * k1] = 0.0;
            py[i * k1 + 1] = 1.0; // first order
            pyOrig[i * k1] = 0.0;
            pyOrig[i * k1 + 1] = 1.0; // first order
        }

        for (size_t j = 0; j < n; j++) {
            x_p[j] = 1;
            x_pOrig[j] = 1;
            tx[j * k1 + 1] = 1;

            _funInner->Forward(1, x_pOrig);
            vector<CGD> dwOrig = _funInner->Reverse(2, pyOrig);
            vector<double> dwInner = modelLib->ReverseTwo(tx, ty, py);
            f2.Forward(1, x_p);
            vector<double> dwOuter = f2.Reverse(2, py);

            x_p[j] = 0;
            x_pOrig[j] = 0;
            tx[j * k1 + 1] = 0;

            // only compare second order information
            // (location of the elements is different then if py.size() == m)
            ASSERT_EQ(dwOrig.size(), n * k1);
            ASSERT_EQ(dwOrig.size(), dwInner.size());
            ASSERT_EQ(dwOrig.size(), dwOuter.size());
            for (size_t j = 0; j < n; j++) {
                ASSERT_TRUE(nearEqual(dwInner[j * k1], dwOrig[j * k1].getValue()));
                ASSERT_TRUE(nearEqual(dwOuter[j * k1], dwOrig[j * k1].getValue()));
            }
        }

        /**
         * Jacobian
         */
        vector<CGD> jacOrig = _funInner->Jacobian(xOrig);
        vector<double> jacOuter = f2.Jacobian(x);
        ASSERT_TRUE(compareValues(jacOuter, jacOrig, epsilonR, epsilonA));

        /**
         * Jacobian sparsity
         */
        const std::vector<bool> jacSparsityOrig = jacobianForwardSparsity<std::vector<bool>, CGD> (*_funInner);
        const std::vector<bool> jacSparsityOuter = jacobianForwardSparsity<std::vector<bool>, double> (f2);

        compareBoolValues(jacSparsityOrig, jacSparsityOuter);

        const std::vector<bool> jacSparsityOrigRev = jacobianReverseSparsity<std::vector<bool>, CGD> (*_funInner);
        const std::vector<bool> jacSparsityOuterRev = jacobianReverseSparsity<std::vector<bool>, double> (f2);

        compareBoolValues(jacSparsityOrigRev, jacSparsityOrig);
        compareBoolValues(jacSparsityOrigRev, jacSparsityOuterRev);

        /**
         * Sparse jacobian
         */
        jacOrig = _funInner->SparseJacobian(xOrig);
        jacOuter = f2.SparseJacobian(x);
        ASSERT_TRUE(compareValues(jacOuter, jacOrig, epsilonR, epsilonA));

        // sparse reverse
        std::vector<size_t> row, col;
        generateSparsityIndexes(jacSparsityOuter, m, n, row, col);

        sparse_jacobian_work workOrig;
        jacOrig.resize(row.size());
        _funInner->SparseJacobianReverse(xOrig, jacSparsityOrig, row, col, jacOrig, workOrig);

        sparse_jacobian_work work2;
        jacOuter.resize(row.size());
        f2.SparseJacobianReverse(x, jacSparsityOuter, row, col, jacOuter, work2);

        ASSERT_TRUE(compareValues(jacOuter, jacOrig, epsilonR, epsilonA));

        /**
         * Hessian
         */
        for (size_t i = 0; i < m; i++) {
            w[i] = 1;
            wOrig[i] = 1;
        }
        vector<CGD> hessOrig = _funInner->Hessian(xOrig, wOrig);
        vector<double> hessOuter = f2.Hessian(x, w);
        ASSERT_TRUE(compareValues(hessOuter, hessOrig, epsilonR, epsilonA));

        /**
         * Sparse Hessian
         */
        hessOrig = _funInner->SparseHessian(xOrig, wOrig);
        hessOuter = f2.SparseHessian(x, w);
        ASSERT_TRUE(compareValues(hessOuter, hessOrig, epsilonR, epsilonA));
    }

    void testAtomicLibAtomicLib(const CppAD::vector<Base>& x,
                                Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        CppAD::vector<Base> xNorm(x.size());
        for (size_t i = 0; i < xNorm.size(); i++)
            xNorm[i] = 1.0;
        CppAD::vector<Base> eqNorm;

        testAtomicLibAtomicLib(x, xNorm, eqNorm, epsilonR, epsilonA);
    }

    /**
     * Test 2 models in 2 dynamic libraries
     */
    void testAtomicLibAtomicLib(const CppAD::vector<Base>& x,
                                const CppAD::vector<Base>& xNorm,
                                const CppAD::vector<Base>& eqNorm,
                                Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        using namespace std;

        prepareAtomicLibAtomicLib(x, xNorm, eqNorm);
        ASSERT_TRUE(_modelLib != nullptr);

        unique_ptr<GenericModel<Base> > modelLibOuter(_dynamicLib2->model(_modelName + "_outer"));
        ASSERT_TRUE(modelLibOuter.get() != nullptr);

        test2LevelAtomicLibModel(_modelLib.get(), modelLibOuter.get(),
                                 x, xNorm, eqNorm, epsilonR, epsilonA);
    }

    /**
     * Test 2 models in the same dynamic library
     */
    void testAtomicLibModelBridge(const CppAD::vector<Base>& x,
                                  const CppAD::vector<Base>& xNorm,
                                  const CppAD::vector<Base>& eqNorm,
                                  Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        using namespace std;

        prepareAtomicLibModelBridge(x, xNorm, eqNorm);

        unique_ptr<GenericModel<Base> > modelLib(_dynamicLib->model(_modelName));
        unique_ptr<GenericModel<Base> > modelLibOuter(_dynamicLib->model(_modelName + "_outer"));

        test2LevelAtomicLibModel(modelLib.get(), modelLibOuter.get(),
                                 x, xNorm, eqNorm, epsilonR, epsilonA);
    }

    /**
     * Test 2 models in the same dynamic library
     */
    void testAtomicLibModelBridgeCustom(const CppAD::vector<Base>& x,
                                        const CppAD::vector<Base>& xNorm,
                                        const CppAD::vector<Base>& eqNorm,
                                        const std::vector<std::set<size_t> >& jacInner,
                                        const std::vector<std::set<size_t> >& hessInner,
                                        const std::vector<std::set<size_t> >& jacOuter,
                                        const std::vector<std::set<size_t> >& hessOuter,
                                        bool createOuterReverse2,
                                        Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        using namespace std;

        prepareAtomicLibModelBridge(x, xNorm, eqNorm,
                                    jacInner, hessInner,
                                    jacOuter, hessOuter,
                                    createOuterReverse2);

        unique_ptr<GenericModel<Base> > modelLib(_dynamicLib->model(_modelName));
        unique_ptr<GenericModel<Base> > modelLibOuter(_dynamicLib->model(_modelName + "_outer"));

        test2LevelAtomicLibModelCustomEls(modelLib.get(), modelLibOuter.get(),
                                          x, xNorm, eqNorm,
                                          jacOuter, hessOuter,
                                          epsilonR, epsilonA);
    }

    /**
     * Test the Jacobian and Hessian sparsity patterns computed directly with
     * the methods of the CGAbstractAtomicFun.
     */
    void testAtomicSparsities(const CppAD::vector<Base>& x) {
        CGAtomicFunBridge<double> atomicfun("innerModel", *_funInner, true);

        //const size_t n = _funInner->Domain();
        const size_t m = _funInner->Range();

        CppAD::vector<CGD> xx(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            xx[i] = x[i];
        }

        CppAD::vector<std::set<size_t>> jacOrig;
        CppAD::vector<std::set<size_t>> jacAtom;

        jacOrig = jacobianForwardSparsitySet<CppAD::vector<std::set<size_t>>, CGD> (*_funInner);
        jacAtom = atomicfun.jacobianForwardSparsitySet(m, xx);
        compareVectorSetValues(jacOrig, jacAtom);

        jacOrig = jacobianReverseSparsitySet<CppAD::vector<std::set<size_t>>, CGD> (*_funInner);
        jacAtom = atomicfun.jacobianReverseSparsitySet(m, xx);
        compareVectorSetValues(jacOrig, jacAtom);

        CppAD::vector<std::set<size_t>> hessOrig;
        CppAD::vector<std::set<size_t>> hessAtom;

        hessOrig = hessianSparsitySet<CppAD::vector<std::set<size_t>>, CGD> (*_funInner);
        hessAtom = atomicfun.hessianSparsitySet(m, xx);
        compareVectorSetValues(hessOrig, hessAtom);
    }

private:

    void test2LevelAtomicLibModel(GenericModel<Base>* modelLib,
                                  GenericModel<Base>* modelLibOuter,
                                  const CppAD::vector<Base>& x,
                                  const CppAD::vector<Base>& xNorm,
                                  const CppAD::vector<Base>& eqNorm,
                                  Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        ASSERT_EQ(x.size(), xNorm.size());

        using namespace CppAD;
        using namespace std;
        using CppAD::vector;

        const size_t n = _funOuter->Domain();
        const size_t m = _funOuter->Range();

        modelLibOuter->addAtomicFunction(modelLib->asAtomic());
        //modelLibOuter->addExternalModel(*modelLib);


        /**
         * Test zero order
         */
        vector<CGD> xOrig(n);
        for (size_t j = 0; j < n; j++)
            xOrig[j] = x[j];

        vector<CGD> yOrig = _funOuter->Forward(0, xOrig);
        vector<double> yOuter = modelLibOuter->ForwardZero(x);

        ASSERT_TRUE(compareValues(yOuter, yOrig, epsilonR, epsilonA));

        /**
         * Test first order forward mode
         */
        size_t k = 1;
        size_t k1 = k + 1;

        vector<double> x_p(n);
        for (size_t j = 0; j < n; j++)
            x_p[j] = 0;

        vector<CGD> x_pOrig(n);
        vector<double> tx(k1 * n);
        for (size_t j = 0; j < n; j++)
            tx[j * k1] = x[j]; // zero order
        for (size_t j = 0; j < n; j++)
            tx[j * k1 + 1] = 0; // first order

        for (size_t j = 0; j < n; j++) {
            x_pOrig[j] = 1;
            tx[j * k1 + 1] = 1;

            vector<CGD> y_pOrig = _funOuter->Forward(1, x_pOrig);
            vector<double> y_pOuter = modelLibOuter->ForwardOne(tx);

            x_pOrig[j] = 0;
            tx[j * k1 + 1] = 0;

            ASSERT_TRUE(compareValues(y_pOuter, y_pOrig, epsilonR, epsilonA));
        }

        /**
         * Test first order reverse mode
         */
        k = 0;
        k1 = k + 1;

        vector<double> w(m);
        for (size_t i = 0; i < m; i++)
            w[i] = 0;
        vector<CGD> wOrig(m);
        tx.resize(k1 * n);
        for (size_t j = 0; j < n; j++)
            tx[j * k1] = x[j]; // zero order
        vector<double> ty(k1 * m);
        for (size_t i = 0; i < m; i++)
            ty[i * k1] = yOuter[i]; // zero order

        for (size_t i = 0; i < m; i++) {
            w[i] = 1;
            wOrig[i] = 1;

            vector<CGD> dwOrig = _funOuter->Reverse(1, wOrig);
            vector<double> dwOuter = modelLibOuter->ReverseOne(tx, ty, w);

            w[i] = 0;
            wOrig[i] = 0;

            ASSERT_TRUE(compareValues(dwOuter, dwOrig, epsilonR, epsilonA));
        }

        /**
         * Test second order reverse mode
         */
        k = 1;
        k1 = k + 1;
        tx.resize(k1 * n);
        ty.resize(k1 * m);
        vector<double> py(k1 * m);
        vector<CGD> pyOrig(k1 * m);
        //wOrig.resize(k1 * m);
        for (size_t j = 0; j < n; j++) {
            tx[j * k1] = x[j]; // zero order
            tx[j * k1 + 1] = 0; // first order
        }
        for (size_t i = 0; i < m; i++) {
            ty[i * k1] = yOuter[i]; // zero order
            py[i * k1] = 0.0;
            py[i * k1 + 1] = 1.0; // first order
            pyOrig[i * k1] = 0.0;
            pyOrig[i * k1 + 1] = 1.0; // first order
        }

        for (size_t j = 0; j < n; j++) {
            x_pOrig[j] = 1;
            tx[j * k1 + 1] = 1;

            _funOuter->Forward(1, x_pOrig);
            vector<CGD> dwOrig = _funOuter->Reverse(2, pyOrig);
            vector<double> dwOuter = modelLibOuter->ReverseTwo(tx, ty, py);

            x_pOrig[j] = 0;
            tx[j * k1 + 1] = 0;

            // only compare second order information
            // (location of the elements is different then if py.size() == m)
            ASSERT_EQ(dwOrig.size(), n * k1);
            ASSERT_EQ(dwOrig.size(), dwOuter.size());
            for (size_t j = 0; j < n; j++) {
                ASSERT_TRUE(nearEqual(dwOuter[j * k1], dwOrig[j * k1].getValue()));
            }
        }

        /**
         * Jacobian
         */
        vector<CGD> jacOrig = _funOuter->Jacobian(xOrig);
        vector<double> jacOuter = modelLibOuter->Jacobian(x);
        ASSERT_TRUE(compareValues(jacOuter, jacOrig, epsilonR, epsilonA));

        /**
         * Jacobian sparsity
         */
        const std::vector<bool> jacSparsityOrig = jacobianSparsity<std::vector<bool>, CGD> (*_funOuter);
        const std::vector<bool> jacSparsityOuter = modelLibOuter->JacobianSparsityBool();

        compareBoolValues(jacSparsityOrig, jacSparsityOuter);

        /**
         * Sparse jacobian
         */
        jacOrig = _funOuter->SparseJacobian(xOrig);
        jacOuter = modelLibOuter->SparseJacobian(x);
        ASSERT_TRUE(compareValues(jacOuter, jacOrig, epsilonR, epsilonA));

        /**
         * Hessian
         */
        for (size_t i = 0; i < m; i++) {
            w[i] = 1;
            wOrig[i] = 1;
        }
        vector<CGD> hessOrig = _funOuter->Hessian(xOrig, wOrig);
        vector<double> hessOuter = modelLibOuter->Hessian(x, w);
        ASSERT_TRUE(compareValues(hessOuter, hessOrig, epsilonR, epsilonA));

        /**
         * Hessian sparsity
         */
        const std::vector<bool> hessSparsityOrig = hessianSparsity<std::vector<bool>, CGD> (*_funOuter);
        const std::vector<bool> hessSparsityOuter = modelLibOuter->HessianSparsityBool();

        compareBoolValues(hessSparsityOrig, hessSparsityOuter);

        /**
         * Sparse Hessian
         */
        hessOrig = _funOuter->SparseHessian(xOrig, wOrig);
        hessOuter = modelLibOuter->SparseHessian(x, w);
        ASSERT_TRUE(compareValues(hessOuter, hessOrig, epsilonR, epsilonA));
    }

    void test2LevelAtomicLibModelCustomEls(GenericModel<Base>* modelLib,
                                           GenericModel<Base>* modelLibOuter,
                                           const CppAD::vector<Base>& x,
                                           const CppAD::vector<Base>& xNorm,
                                           const CppAD::vector<Base>& eqNorm,
                                           const std::vector<std::set<size_t> >& jacOuterSpar,
                                           const std::vector<std::set<size_t> >& hessOuterSpar,
                                           Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        ASSERT_EQ(x.size(), xNorm.size());

        using namespace CppAD;
        using namespace std;
        using CppAD::vector;

        const size_t n = _funOuter->Domain();
        const size_t m = _funOuter->Range();

        modelLibOuter->addAtomicFunction(modelLib->asAtomic());

        /**
         * Test zero order
         */
        vector<CGD> xOrig(n);
        for (size_t j = 0; j < n; j++)
            xOrig[j] = x[j];

        vector<CGD> yOrig = _funOuter->Forward(0, xOrig);
        vector<double> yOuter = modelLibOuter->ForwardZero(x);

        ASSERT_TRUE(compareValues(yOuter, yOrig, epsilonR, epsilonA));

        /**
         * Jacobian sparsity
         */
        const std::vector<bool> jacSparsityOrig = jacobianSparsity<std::vector<bool>, CGD> (*_funOuter);

        /**
         * Sparse jacobian
         */
        CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
        std::vector<size_t> jacOuterRows, jacOuterCols;
        generateSparsityIndexes(jacOuterSpar, jacOuterRows, jacOuterCols);
        vector<CGD> jacOrig(jacOuterCols.size());
        _funOuter->SparseJacobianReverse(xOrig, jacSparsityOrig,
                                     jacOuterRows, jacOuterCols, jacOrig, work);

        std::vector<double> jacOuter(jacOuterRows.size());
        size_t const* rows, *cols;
        modelLibOuter->SparseJacobian(x, jacOuter, &rows, &cols);

        ASSERT_TRUE(compareValues(jacOuter, jacOrig, epsilonR, epsilonA));

        /**
         * Hessian sparsity
         */
        const std::vector<bool> hessSparsityOrig = hessianSparsity<std::vector<bool>, CGD> (*_funOuter);
        /**
        printSparsityPattern(hessSparsityOrig, "original", n, n);

        std::vector<size_t> hessOuterRows2, hessOuterCols2;
        modelLibOuter->HessianSparsity(hessOuterRows2, hessOuterCols2);
        printSparsityPattern(hessOuterRows2, hessOuterCols2, "outer", n);
        */
        
        /**
         * Sparse Hessian
         */
        vector<CGD> wOrig(m);
        vector<double> w(m);
        for (size_t i = 0; i < m; i++) {
            wOrig[i] = 1;
            w[i] = 1;
        }
        CppAD::sparse_hessian_work hessWork;
        std::vector<size_t> hessOuterRows, hessOuterCols;
        generateSparsityIndexes(hessOuterSpar, hessOuterRows, hessOuterCols);
        vector<CGD> hessOrig(hessOuterRows.size());
        _funOuter->SparseHessian(xOrig, wOrig, hessSparsityOrig,
                                 hessOuterRows, hessOuterCols, hessOrig, hessWork);

        vector<double> hessOuter(hessOuterRows.size());
        modelLibOuter->SparseHessian(x, w,
                                     hessOuter, &rows, &cols);

        ASSERT_TRUE(compareValues(hessOuter, hessOrig, epsilonR, epsilonA));
    }

    void prepareAtomicLib(const CppAD::vector<Base>& x,
                          const CppAD::vector<Base>& xNorm,
                          const CppAD::vector<Base>& eqNorm) {
        const size_t n = x.size();

        // independent variables
        std::vector<ADCGD> u(n);
        for (size_t j = 0; j < n; j++)
            u[j] = x[j];

        CppAD::Independent(u);
        for (size_t j = 0; j < n; j++)
            u[j] *= xNorm[j];

        /**
         * create the CppAD tape as usual
         */
        std::vector<ADCGD> Z = model(u);

        if (eqNorm.size() > 0) {
            ASSERT_EQ(Z.size(), eqNorm.size());
            for (size_t i = 0; i < Z.size(); i++)
                Z[i] /= eqNorm[i];
        }

        // create f: U -> Z and vectors used for derivative calculations
        _funInner.reset(new ADFun<CGD>());
        _funInner->Dependent(Z);

        /**
         * Create the dynamic library
         * (generate and compile source code)
         */
        ModelCSourceGen<double> compHelp(*_funInner, _modelName);

        compHelp.setCreateForwardZero(true);
        compHelp.setCreateForwardOne(true);
        compHelp.setCreateReverseOne(true);
        compHelp.setCreateReverseTwo(true);

        ModelLibraryCSourceGen<double> compDynHelp(compHelp);
        compDynHelp.setVerbose(this->verbose_);

        SaveFilesModelLibraryProcessor<double>::saveLibrarySourcesTo(compDynHelp, "sources_atomiclib_" + _modelName);

        DynamicModelLibraryProcessor<double> p(compDynHelp);

        GccCompiler<double> compiler;
        prepareTestCompilerFlags(compiler);
        _dynamicLib.reset(p.createDynamicLibrary(compiler));
    }

    virtual void prepareAtomicLibAtomicLib(const CppAD::vector<Base>& x,
                                           const CppAD::vector<Base>& xNorm,
                                           const CppAD::vector<Base>& eqNorm) {
        using CppAD::vector;

        const size_t n = x.size();

        // independent variables
        std::vector<ADCGD> ax(n);
        for (size_t j = 0; j < n; j++)
            ax[j] = x[j];

        CppAD::Independent(ax);
        for (size_t j = 0; j < n; j++)
            ax[j] *= xNorm[j];

        /**
         * create the CppAD tape as usual
         */
        std::vector<ADCGD> ay = model(ax);

        if (eqNorm.size() > 0) {
            ASSERT_EQ(ay.size(), eqNorm.size());
            for (size_t i = 0; i < ay.size(); i++)
                ay[i] /= eqNorm[i];
        }

        // create f: U -> Z and vectors used for derivative calculations
        _funInner.reset(new ADFun<CGD>());
        _funInner->Dependent(ay);

        /**
         * Create the dynamic library model
         */
        ModelCSourceGen<double> cSourceInner(*_funInner, _modelName);
        cSourceInner.setCreateForwardZero(true);
        cSourceInner.setCreateForwardOne(true);
        cSourceInner.setCreateReverseOne(true);
        cSourceInner.setCreateReverseTwo(true);

        /**
         * Create the dynamic library
         * (generate and compile source code)
         */
        GccCompiler<double> compiler1;
        prepareTestCompilerFlags(compiler1);

        ModelLibraryCSourceGen<double> compDynHelp(cSourceInner);
        compDynHelp.setVerbose(this->verbose_);

        SaveFilesModelLibraryProcessor<double>::saveLibrarySourcesTo(compDynHelp, "sources_atomiclibatomiclib_" + _modelName);

        DynamicModelLibraryProcessor<double> p(compDynHelp, "innerModel");
        _dynamicLib.reset(p.createDynamicLibrary(compiler1));
        _modelLib.reset(_dynamicLib->model(_modelName));

        /**
         * Second compiled model
         */
        // independent variables
        std::vector<ADCGD> ax2(n);
        for (size_t j = 0; j < n; j++)
            ax2[j] = x[j];

        CppAD::Independent(ax2);

        CGAtomicGenericModel<Base>& innerAtomicFun = _modelLib->asAtomic();
        CGAtomicFun<Base> cgInnerAtomicFun(innerAtomicFun, std::vector<double>(_modelLib->Domain()), true); // required for taping

        std::vector<ADCGD> ZZ(ay.size());
        cgInnerAtomicFun(ax2, ZZ);
        std::vector<ADCGD> Z2 = modelOuter(ZZ);

        // create f: U2 -> Z2
        ADFun<CGD> fun2;
        fun2.Dependent(Z2);

        ModelCSourceGen<double> cSourceOuter(fun2, _modelName + "_outer");
        cSourceOuter.setCreateForwardZero(true);
        cSourceOuter.setCreateForwardOne(true);
        cSourceOuter.setCreateReverseOne(true);
        cSourceOuter.setCreateReverseTwo(true);
        cSourceOuter.setCreateJacobian(true);
        cSourceOuter.setCreateHessian(true);
        cSourceOuter.setCreateSparseJacobian(true);
        cSourceOuter.setCreateSparseHessian(true);
        //cSourceOuter.setMaxAssignmentsPerFunc(20);

        /**
         * Create the dynamic library
         * (generate and compile source code)
         */
        ModelLibraryCSourceGen<double> compDynHelp2(cSourceOuter);
        compDynHelp2.setVerbose(this->verbose_);

        SaveFilesModelLibraryProcessor<double>::saveLibrarySourcesTo(compDynHelp2, "sources_atomiclibatomiclib_" + _modelName);

        DynamicModelLibraryProcessor<double> p2(compDynHelp2, "outerModel");
        GccCompiler<double> compiler2;
        prepareTestCompilerFlags(compiler2);
        _dynamicLib2.reset(p2.createDynamicLibrary(compiler2));

        /**
         * tape the model without atomics
         */
        tapeOuterModel(x, xNorm, eqNorm);
    }

    virtual void prepareAtomicLibModelBridge(const CppAD::vector<Base>& x,
                                             const CppAD::vector<Base>& xNorm,
                                             const CppAD::vector<Base>& eqNorm) {
        std::vector<std::set<size_t> > jacInner, hessInner;
        std::vector<std::set<size_t> > jacOuter, hessOuter;
        prepareAtomicLibModelBridge(x, xNorm, eqNorm,
                                    jacInner, hessInner,
                                    jacOuter, hessOuter,
                                    true);
    }

    virtual void prepareAtomicLibModelBridge(const CppAD::vector<Base>& x,
                                             const CppAD::vector<Base>& xNorm,
                                             const CppAD::vector<Base>& eqNorm,
                                             const std::vector<std::set<size_t> >& jacInner,
                                             const std::vector<std::set<size_t> >& hessInner,
                                             const std::vector<std::set<size_t> >& jacOuter,
                                             const std::vector<std::set<size_t> >& hessOuter,
                                             bool createOuterReverse2) {
        const size_t n = x.size();

        // independent variables
        std::vector<ADCGD> u(n);
        for (size_t j = 0; j < n; j++)
            u[j] = x[j];

        CppAD::Independent(u);
        for (size_t j = 0; j < n; j++)
            u[j] *= xNorm[j];

        /**
         * create the CppAD tape as usual
         */
        std::vector<ADCGD> Z = model(u);

        if (eqNorm.size() > 0) {
            ASSERT_EQ(Z.size(), eqNorm.size());
            for (size_t i = 0; i < Z.size(); i++)
                Z[i] /= eqNorm[i];
        }

        // create f: U -> Z and vectors used for derivative calculations
        _funInner.reset(new ADFun<CGD>());
        _funInner->Dependent(Z);

        /**
         * Create the dynamic library model
         */
        ModelCSourceGen<double> cSourceInner(*_funInner, _modelName);
        if (jacInner.size() > 0) {
            cSourceInner.setCustomSparseJacobianElements(jacInner);
        }
        if (hessInner.size() > 0) {
            cSourceInner.setCustomSparseHessianElements(hessInner);
        }

        /**
         * Second compiled model
         */
        // independent variables
        std::vector<ADCGD> u2(n);
        for (size_t j = 0; j < n; j++)
            u2[j] = x[j];

        CppAD::Independent(u2);

        CGAtomicFunBridge<double> atomicfun(_modelName, *_funInner, true);
        if (jacInner.size() > 0) {
            atomicfun.setCustomSparseJacobianElements(jacInner);
        }
        if (hessInner.size() > 0) {
            atomicfun.setCustomSparseHessianElements(hessInner);
        }

        std::vector<ADCGD> ZZ(Z.size());
        atomicfun(u2, ZZ);
        std::vector<ADCGD> Z2 = modelOuter(ZZ);

        // create f: U2 -> Z2
        ADFun<CGD> fun2;
        fun2.Dependent(Z2);

        ModelCSourceGen<double> cSourceOuter(fun2, _modelName + "_outer");

        cSourceInner.setCreateForwardZero(true);
        cSourceInner.setCreateForwardOne(true);
        cSourceInner.setCreateReverseOne(true);
        cSourceInner.setCreateReverseTwo(true);
        //cSourceInner.setCreateSparseHessian(true); //not really required

        cSourceOuter.setCreateForwardZero(true);
        cSourceOuter.setCreateForwardOne(true);
        cSourceOuter.setCreateReverseOne(true);
        cSourceOuter.setCreateReverseTwo(createOuterReverse2);
        cSourceOuter.setCreateJacobian(true);
        cSourceOuter.setCreateHessian(true);
        cSourceOuter.setCreateSparseJacobian(true);
        cSourceOuter.setCreateSparseHessian(true);
        if (jacOuter.size() > 0) {
            cSourceOuter.setCustomSparseJacobianElements(jacOuter);
        }
        if (hessOuter.size() > 0) {
            cSourceOuter.setCustomSparseHessianElements(hessOuter);
        }

        /**
         * Create the dynamic library
         * (generate and compile source code)
         */
        ModelLibraryCSourceGen<double> compDynHelp(cSourceInner, cSourceOuter);
        compDynHelp.setVerbose(this->verbose_);

        std::string folder = std::string("sources_atomiclibmodelbridge_") + (createOuterReverse2 ? "rev2_" : "dir_") + _modelName;
        SaveFilesModelLibraryProcessor<double>::saveLibrarySourcesTo(compDynHelp, folder);

        DynamicModelLibraryProcessor<double> p(compDynHelp);

        GccCompiler<double> compiler;
        prepareTestCompilerFlags(compiler);
        _dynamicLib.reset(p.createDynamicLibrary(compiler));

        /**
         * tape the model without atomics
         */
        tapeOuterModel(x, xNorm, eqNorm);
    }

    void tapeOuterModel(const CppAD::vector<Base>& x,
                        const CppAD::vector<Base>& xNorm,
                        const CppAD::vector<Base>& eqNorm) {
        const size_t n = x.size();

        // independent variables
        std::vector<ADCGD> u(n);
        for (size_t j = 0; j < n; j++)
            u[j] = x[j];

        CppAD::Independent(u);
        for (size_t j = 0; j < n; j++)
            u[j] *= xNorm[j];

        /**
         * create the CppAD tape as usual
         */
        std::vector<ADCGD> yInner = model(u);
        if (eqNorm.size() > 0) {
            ASSERT_EQ(yInner.size(), eqNorm.size());
            for (size_t i = 0; i < yInner.size(); i++)
                yInner[i] /= eqNorm[i];
        }

        std::vector<ADCGD> yOuter = modelOuter(yInner);

        // create f: U -> Z2
        _funOuter.reset(new ADFun<CGD>());
        _funOuter->Dependent(yOuter);
    }

};

} // END cg namespace
} // END CppAD namespace

#endif
