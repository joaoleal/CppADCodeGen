#ifndef CPPAD_CG_TEST_CPPADCGDYNAMICATOMICTEST_INCLUDED
#define CPPAD_CG_TEST_CPPADCGDYNAMICATOMICTEST_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
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

#include <memory>

#include "CppADCGModelTest.hpp"
#include "gccCompilerFlags.hpp"

namespace CppAD {
namespace cg {

class CppADCGDynamicAtomicTest : public CppADCGModelTest {
public:
    using Super = CppADCGTest;
    using Base = Super::Base;
    using CGD = Super::CGD;
    using ADCGD = Super::ADCGD;
protected:
    const std::string _modelName;
    std::unique_ptr<ADFun<CGD>> _funInner; // inner model tape
    std::unique_ptr<ADFun<CGD>> _funOuterNoAtom; // outer model tape with no atomics
    std::unique_ptr<ADFun<CGD>> _funOuterAtom; // outer model tape with atomics
    std::unique_ptr<ADFun<Base>> _funBaseOuterAtom; // outer model tape with atomics
    std::unique_ptr<DynamicLib<Base>> _dynamicLib; // library for the inner model (in may also have the outer model)
    std::unique_ptr<DynamicLib<Base>> _dynamicLibOuter; // library for the outer model
    std::unique_ptr<GenericModel<Base>> _modelLib;
    std::unique_ptr<CGAtomicFun<Base>> _atomFun;
public:

    explicit CppADCGDynamicAtomicTest(std::string modelName,
                                      bool verbose = false,
                                      bool printValues = false) :
            CppADCGModelTest(verbose, printValues),
            _modelName(std::move(modelName)) {
        //this->verbose_ = true;
    }

    virtual std::vector<ADCGD> model(const std::vector<ADCGD>& x,
                                     const std::vector<ADCGD>& par) = 0;

    virtual std::vector<ADCGD> modelOuter(const std::vector<ADCGD>& yInner,
                                          const std::vector<ADCGD>& par) {
        std::vector<ADCGD> y(yInner.size() - 1);

        for (size_t i = 0; i < yInner.size() - 1; i++) {
            y[i] = 2 * yInner[i];
        }
        y[y.size() - 1] += yInner[yInner.size() - 1];

        return y;
    }

    void TearDown() override {
        // need to clear memory here because CppAD checks if there is any memory still in use
        _dynamicLib.reset();
        _dynamicLibOuter.reset();
        _modelLib.reset();
        _funInner.reset();
        _funOuterNoAtom.reset();
        _funOuterAtom.reset();
        _funBaseOuterAtom.reset();
        _atomFun.reset();
    }

    ~CppADCGDynamicAtomicTest() override = default;

    inline static std::vector<Base> makeValues(size_t n) {
        std::vector<Base> x(n);

        for (size_t j = 0; j < n; j++)
            x[j] = j + 2;

        return x;
    }

    /**
     * Tests one compiled model used as an atomic function by an ADFun.
     * The outer model method (modelOuter) is NOT used here!
     * It uses z = g(x) = f(x).
     */
    void testADFunAtomicLibSimple(const std::vector<Base>& x,
                                  const std::vector<Base>& par = {},
                                  const std::vector<Base>& xNorm = {},
                                  const std::vector<Base>& eqNorm = {},
                                  Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        ASSERT_TRUE(xNorm.empty() || x.size() == xNorm.size());

        using namespace CppAD;
        using namespace std;
        using CppAD::vector;

        tapeInnerModel(x, par, xNorm, eqNorm);

        unique_ptr<GenericModel<Base> > modelLib = compileInnerModel("");

        const size_t m = _funInner->Range();

        auto ax = makeADVector(x);
        auto ap = makeADVector(par);

        // declare independent variables and start tape recording
        size_t abort_op_index = 0;
        bool record_compare = true;
        CppAD::Independent(ax, abort_op_index, record_compare, ap);

        std::vector<AD<double> > ay(m);

        // call user function and store CGAtomicLibModel(x) in au[0]
        CGAtomicGenericModel<double>& atomicfun = modelLib->asAtomic();

        atomicfun(ax, ay);

        // create fWrapAtom: x -> y and stop tape recording
        ADFun<double> fWrapAtom(ax, ay);

        /**
         * Test zero order
         */
         testForwardZeroResults(*modelLib, *_funInner, &fWrapAtom, x, par, epsilonR, epsilonA);

        /**
         * Test first order forward mode
         */
        testForwardOneResults(*modelLib, *_funInner, &fWrapAtom, x, par, epsilonR, epsilonA);

        /**
         * Test first order reverse mode
         */
        testReverseOneResults(*modelLib, *_funInner, &fWrapAtom, x, par, epsilonR, epsilonA);

        /**
         * Test second order reverse mode
         */
        testReverseTwoResults(*modelLib, *_funInner, &fWrapAtom, x, par, epsilonR, epsilonA);

        /**
         * Dense Jacobian
         */
        testDenseJacResults(*modelLib, *_funInner, x, par, epsilonR, epsilonA);

        /**
         * Jacobian sparsity
         */
        testJacobianSparsity(*_funInner, fWrapAtom);

        /**
         * Sparse Jacobian
         */
        size_t n_tests = _dynamicLib->getThreadNumber() > 1 ? 2 : 1;
        testSparseJacobianResults(n_tests, *modelLib, *_funInner, &fWrapAtom, x, par, false, epsilonR, epsilonA);

        /**
         * Dense Hessian
         */
        testDenseHessianResults(*modelLib, *_funInner, x, par, epsilonR, epsilonA);

        /**
         * Sparse Hessian
         */
        testSparseHessianResults(n_tests, *modelLib, *_funInner, &fWrapAtom, x, par, false, epsilonR, epsilonA);
    }

    /**
     * Tests one compiled model used as an atomic function by an ADFun.
     * The outer model method (modelOuter) is used!
     * It uses z = g(f(x)).
     */
    void testADFunAtomicLib(const std::vector<Base>& x,
                            const std::vector<Base>& par,
                            const std::vector<Base>& xNorm = {},
                            const std::vector<Base>& eqNorm = {},
                            Base epsilonR = 1e-14,
                            Base epsilonA = 1e-14) {
        prepareAtomicLibAtomicLib(x, par, xNorm, eqNorm);
        ASSERT_TRUE(_modelLib != nullptr);
        ASSERT_TRUE(_funOuterNoAtom != nullptr);

        auto atomic = [&](const std::vector<ADCGD>& ax,
                          std::vector<ADCGD>& ay){
            assert(_modelLib != nullptr);
            assert(_modelLib->getName() ==  _modelName);
            _atomFun = std::make_unique<CGAtomicFun<double>>(_modelLib->asAtomic(), x, par.size());

            (*_atomFun)(ax, ay); // atomic function
        };

        tapeOuterModelWithAtomic(atomic, x, par, xNorm, eqNorm);
        ASSERT_TRUE(_funOuterAtom != nullptr);

        testAtomicLibModelInCppAD(*_funOuterNoAtom, *_funOuterAtom, x, par, epsilonR, epsilonA);
    }

    /**
     * Test 2 models in 2 dynamic libraries
     */
    void testAtomicLibAtomicLib(const std::vector<Base>& x,
                                const std::vector<Base>& par,
                                const std::vector<Base>& xNorm = {},
                                const std::vector<Base>& eqNorm = {},
                                Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        using namespace std;

        prepareAtomicLibAtomicLib(x, par, xNorm, eqNorm);
        ASSERT_TRUE(_modelLib != nullptr);

        unique_ptr<GenericModel<Base> > modelLibOuter = _dynamicLibOuter->model(_modelName + "_outer");
        ASSERT_TRUE(modelLibOuter != nullptr);

        test2LevelAtomicLibModel(_modelLib.get(), modelLibOuter.get(),
                                 x, par, xNorm, eqNorm, epsilonR, epsilonA);
    }

    /**
     * Test 2 models in the same dynamic library
     */
    void testAtomicLibModelBridge(const std::vector<Base>& x,
                                  const std::vector<Base>& par,
                                  const std::vector<Base>& xNorm,
                                  const std::vector<Base>& eqNorm,
                                  Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        using namespace std;

        prepareAtomicLibModelBridge(x, par, xNorm, eqNorm);

        unique_ptr<GenericModel<Base> > modelLib = _dynamicLib->model(_modelName);
        unique_ptr<GenericModel<Base> > modelLibOuter = _dynamicLib->model(_modelName + "_outer");

        test2LevelAtomicLibModel(modelLib.get(), modelLibOuter.get(),
                                 x, par, xNorm, eqNorm, epsilonR, epsilonA);
    }

    /**
     * Test 2 models in the same dynamic library
     */
    void testAtomicLibModelBridgeCustom(const std::vector<Base>& x,
                                        const std::vector<Base>& par,
                                        const std::vector<Base>& xNorm,
                                        const std::vector<Base>& eqNorm,
                                        const std::vector<std::set<size_t> >& jacInner,
                                        const std::vector<std::set<size_t> >& hessInner,
                                        const std::vector<std::set<size_t> >& jacOuter,
                                        const std::vector<std::set<size_t> >& hessOuter,
                                        bool createOuterReverse2,
                                        Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        using namespace std;

        prepareAtomicLibModelBridge(x, par, xNorm, eqNorm,
                                    jacInner, hessInner,
                                    jacOuter, hessOuter,
                                    createOuterReverse2);

        unique_ptr<GenericModel<Base> > modelLib = _dynamicLib->model(_modelName);
        unique_ptr<GenericModel<Base> > modelLibOuter = _dynamicLib->model(_modelName + "_outer");

        test2LevelAtomicLibModelCustomEls(modelLib.get(), modelLibOuter.get(),
                                          x, par, xNorm, eqNorm,
                                          jacOuter, hessOuter,
                                          epsilonR, epsilonA);
    }

    /**
     * Test the Jacobian and Hessian sparsity patterns computed directly with
     * the methods of the CGAbstractAtomicFun.
     */
    void testAtomicSparsities(const std::vector<Base>& x,
                              const std::vector<Base>& par) {
        CGAtomicFunBridge<double> atomicfun("innerModel", *_funInner, true);

        //const size_t n = _funInner->Domain();
        //const size_t m = _funInner->Range();

        CppAD::vector<CGD> pp(x.size() + par.size());
        std::copy(x.data(), x.data() + x.size(), pp.data());
        std::copy(par.data(), par.data() + par.size(), pp.data() + x.size());

        CppAD::vector<std::set<size_t>> jacOrig;
        CppAD::vector<std::set<size_t>> jacAtom;

        jacOrig = jacobianForwardSparsitySet<CppAD::vector<std::set<size_t>>, CGD> (*_funInner);
        generateSparsitySet(atomicfun.jacobianSparsity(pp), jacAtom);
        compareVectorSetValues(jacOrig, jacAtom);

        CppAD::vector<std::set<size_t>> hessOrig;
        CppAD::vector<std::set<size_t>> hessAtom;

        hessOrig = hessianSparsitySet<CppAD::vector<std::set<size_t>>, CGD> (*_funInner);
        generateSparsitySet(atomicfun.hessianSparsity(pp), hessAtom);
        compareVectorSetValues(hessOrig, hessAtom);
    }

private:

    void test2LevelAtomicLibModel(GenericModel<Base>* modelLib,
                                  GenericModel<Base>* modelLibOuter,
                                  const std::vector<Base>& x,
                                  const std::vector<Base>& par,
                                  const std::vector<Base>& xNorm,
                                  const std::vector<Base>& eqNorm,
                                  Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        ASSERT_TRUE(xNorm.empty() || x.size() == xNorm.size());

        using namespace CppAD;
        using namespace std;
        using CppAD::vector;

        std::cout << std::endl << "2LevelAtomicLib" << std::endl;

        const size_t n = _funOuterNoAtom->Domain();
        const size_t m = _funOuterNoAtom->Range();

        modelLibOuter->addAtomicFunction(modelLib->asAtomic());

        /**
         * Test zero order
         */
        std::vector<CGD> xOrig(n);
        for (size_t j = 0; j < n; j++)
            xOrig[j] = x[j];

        auto yOrig = _funOuterNoAtom->Forward(0, xOrig);
        auto yOuter = modelLibOuter->ForwardZero(x);

        ASSERT_TRUE(compareValues<double>(yOuter, yOrig, epsilonR, epsilonA));

        /**
         * Test first order forward mode
         */
        size_t k = 1;
        size_t k1 = k + 1;

        vector<double> x_p(n);
        for (size_t j = 0; j < n; j++)
            x_p[j] = 0;

        vector<CGD> x_pOrig(n);
        std::vector<double> tx(k1 * n);
        for (size_t j = 0; j < n; j++)
            tx[j * k1] = x[j]; // zero order
        for (size_t j = 0; j < n; j++)
            tx[j * k1 + 1] = 0; // first order

        for (size_t j = 0; j < n; j++) {
            x_pOrig[j] = 1;
            tx[j * k1 + 1] = 1;

            auto y_pOrig = _funOuterNoAtom->Forward(1, x_pOrig);
            auto y_pOuter = modelLibOuter->ForwardOne(tx, par);

            x_pOrig[j] = 0;
            tx[j * k1 + 1] = 0;

            ASSERT_TRUE(compareValues<double>(y_pOuter, y_pOrig, epsilonR, epsilonA));
        }

        /**
         * Test first order reverse mode
         */
        k = 0;
        k1 = k + 1;

        std::vector<double> w(m);
        for (size_t i = 0; i < m; i++)
            w[i] = 0;
        std::vector<CGD> wOrig(m);
        tx.resize(k1 * n);
        for (size_t j = 0; j < n; j++)
            tx[j * k1] = x[j]; // zero order
        std::vector<double> ty(k1 * m);
        for (size_t i = 0; i < m; i++)
            ty[i * k1] = yOuter[i]; // zero order

        for (size_t i = 0; i < m; i++) {
            w[i] = 1;
            wOrig[i] = 1;

            auto dwOrig = _funOuterNoAtom->Reverse(1, wOrig);
            auto dwOuter = modelLibOuter->ReverseOne(tx, par, ty, w);

            w[i] = 0;
            wOrig[i] = 0;

            ASSERT_TRUE(compareValues<double>(dwOuter, dwOrig, epsilonR, epsilonA));
        }

        /**
         * Test second order reverse mode
         */
        k = 1;
        k1 = k + 1;
        tx.resize(k1 * n);
        ty.resize(k1 * m);
        std::vector<double> py(k1 * m);
        std::vector<CGD> pyOrig(k1 * m);
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

            _funOuterNoAtom->Forward(1, x_pOrig);
            auto dwOrig = _funOuterNoAtom->Reverse(2, pyOrig);
            auto dwOuter = modelLibOuter->ReverseTwo(tx, par, ty, py);

            x_pOrig[j] = 0;
            tx[j * k1 + 1] = 0;

            // only compare second order information
            // (location of the elements is different then if py.size() == m)
            ASSERT_EQ(dwOrig.size(), n * k1);
            ASSERT_EQ(dwOrig.size(), dwOuter.size());
            for (size_t j2 = 0; j2 < n; j2++) {
                ASSERT_TRUE(nearEqual(dwOuter[j2 * k1], dwOrig[j2 * k1].getValue()));
            }
        }

        /**
         * Jacobian
         */
        auto jacOrig = _funOuterNoAtom->Jacobian(xOrig);
        auto jacOuter = modelLibOuter->Jacobian(x, par);
        ASSERT_TRUE(compareValues<double>(jacOuter, jacOrig, epsilonR, epsilonA));

        /**
         * Jacobian sparsity
         */
        const std::vector<bool> jacSparsityOrig = jacobianSparsityBool<std::vector<bool>, CGD> (*_funOuterNoAtom);
        const std::vector<bool> jacSparsityOuter = modelLibOuter->JacobianSparsityBool();

        compareBoolValues(jacSparsityOrig, jacSparsityOuter);

        /**
         * Sparse jacobian
         */
        jacOrig = _funOuterNoAtom->SparseJacobian(xOrig);
        jacOuter = modelLibOuter->SparseJacobian(x, par);
        ASSERT_TRUE(compareValues<double>(jacOuter, jacOrig, epsilonR, epsilonA));

        /**
         * Hessian
         */
        for (size_t i = 0; i < m; i++) {
            w[i] = 1;
            wOrig[i] = 1;
        }
        auto hessOrig = _funOuterNoAtom->Hessian(xOrig, wOrig);
        auto hessOuter = modelLibOuter->Hessian(x, par, w);
        ASSERT_TRUE(compareValues<double>(hessOuter, hessOrig, epsilonR, epsilonA));

        /**
         * Hessian sparsity
         */
        const std::vector<bool> hessSparsityOrig = hessianSparsityBool<std::vector<bool>, CGD> (*_funOuterNoAtom);
        const std::vector<bool> hessSparsityOuter = modelLibOuter->HessianSparsityBool();

        compareBoolValues(hessSparsityOrig, hessSparsityOuter);

        /**
         * Sparse Hessian
         */
        hessOrig = _funOuterNoAtom->SparseHessian(xOrig, wOrig);
        hessOuter = modelLibOuter->SparseHessian(x, par, w);
        ASSERT_TRUE(compareValues<double>(hessOuter, hessOrig, epsilonR, epsilonA));
    }

    static void testAtomicLibModelInCppAD(ADFun<CGD>& funOuter,
                                          ADFun<CGD>& funOuterAtom, // with an atomic function
                                          const std::vector<Base>& xx,
                                          const std::vector<Base>& par,
                                          Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        using namespace CppAD;
        using namespace std;
        using CppAD::vector;

        std::cout << std::endl << "AtomicLibModelInCppAD" << std::endl;

        const size_t n = funOuter.Domain();
        const size_t m = funOuter.Range();

        /**
         * Test zero order
         */
        auto x = makeCGVector(xx);

        auto y = funOuter.Forward(0, x);
        auto yAtom = funOuterAtom.Forward(0, x); // with an atomic function

        ASSERT_TRUE(compareValues<double>(y, yAtom, epsilonR, epsilonA));

        /**
         * Test first order forward mode
         */

        vector<CGD> x_p(n);

        for (size_t j = 0; j < n; j++) {
            x_p[j] = 1;

            vector<CGD> y_p = funOuter.Forward(1, x_p);
            vector<CGD> y_pAtom = funOuterAtom.Forward(1, x_p); // with an atomic function

            x_p[j] = 0;

            ASSERT_TRUE(compareValues<double>(y_p, y_pAtom, epsilonR, epsilonA));
        }

        /**
         * Test first order reverse mode
         */
        std::vector<CGD> w(m);
        for (size_t i = 0; i < m; i++) {
            w[i] = 1;

            auto dw = funOuter.Reverse(1, w);
            auto dwAtom = funOuterAtom.Reverse(1, w); // with an atomic function

            w[i] = 0;

            ASSERT_TRUE(compareValues<double>(dw, dwAtom, epsilonR, epsilonA));
        }

        /**
         * Test second order reverse mode
         */
        size_t k = 1;
        size_t k1 = k + 1;
        vector<CGD> py(m); // not (k1 * m)
        for (size_t i = 0; i < m; i++) {
            //py[i * k1] = 0.0;
            //py[i * k1 + 1] = 1.0; // first order
            py[i] = 1.0; // first order
        }

        for (size_t j = 0; j < n; j++) {
            x_p[j] = 1;

            funOuter.Forward(1, x_p);
            funOuterAtom.Forward(1, x_p);  // with an atomic function
            vector<CGD> dw = funOuter.Reverse(2, py);
            vector<CGD> dwAtom = funOuterAtom.Reverse(2, py); // with an atomic function

            x_p[j] = 0;

            // only compare second order information
            ASSERT_EQ(dw.size(), n * k1);
            ASSERT_EQ(dw.size(), dwAtom.size());
            for (size_t j2 = 0; j2 < n; j2++) {
                ASSERT_TRUE(nearEqual(dw[j2 * k1 + 1].getValue(), dwAtom[j2 * k1 + 1].getValue()));
            }
        }

        /**
         * Jacobian
         */
        auto jac = funOuter.Jacobian(x);
        auto jacAtom = funOuterAtom.Jacobian(x); // with an atomic function
        ASSERT_TRUE(compareValues<double>(jac, jacAtom, epsilonR, epsilonA));

        /**
         * Jacobian sparsity
         */
        const std::vector<bool> jacSparsity = jacobianSparsityBool<std::vector<bool>, CGD> (funOuter);
        const std::vector<bool> jacSparsityAtom = jacobianSparsityBool<std::vector<bool>, CGD> (funOuterAtom);

        compareBoolValues(jacSparsity, jacSparsityAtom);

        /**
         * Sparse jacobian
         */
        jac = funOuter.SparseJacobian(x);
        jacAtom = funOuterAtom.SparseJacobian(x); // with an atomic function
        ASSERT_TRUE(compareValues<double>(jac, jacAtom, epsilonR, epsilonA));

        /**
         * Hessian
         */
        for (size_t i = 0; i < m; i++) {
            w[i] = 1;
        }
        auto hess = funOuter.Hessian(x, w);
        auto hessAtom = funOuterAtom.Hessian(x, w); // with an atomic function
        ASSERT_TRUE(compareValues<double>(hess, hessAtom, epsilonR, epsilonA));

        /**
         * Hessian sparsity
         */
        const std::vector<bool> hessSparsity = hessianSparsityBool<std::vector<bool>, CGD> (funOuter);
        const std::vector<bool> hessSparsityAtom = hessianSparsityBool<std::vector<bool>, CGD> (funOuterAtom);

        compareBoolValues(hessSparsity, hessSparsityAtom);

        /**
         * Sparse Hessian
         */
        hess = funOuter.SparseHessian(x, w);
        hessAtom = funOuterAtom.SparseHessian(x, w); // with an atomic function
        ASSERT_TRUE(compareValues<double>(hess, hessAtom, epsilonR, epsilonA));

        //CodeHandler<Base> handler;
        //handler.makeVariables(x);
        //handler.makeVariables(w);
        //hessAtom = funOuterAtom.SparseHessian(x, w); // with an atomic function
        //
        //LanguageC<Base> lang("double");
        //LangCDefaultVariableNameGenerator<Base> nameGen("hess");
        //LangCDefaultHessianVarNameGenerator<Base> nameGenHess(&nameGen, n);
        //handler.generateCode(std::cout, lang, hessAtom, nameGenHess);
        //std::cout << std::endl;
    }

    void test2LevelAtomicLibModelCustomEls(GenericModel<Base>* modelLib,
                                           GenericModel<Base>* modelLibOuter,
                                           const std::vector<Base>& x,
                                           const std::vector<Base>& par,
                                           const std::vector<Base>& xNorm,
                                           const std::vector<Base>& eqNorm,
                                           const std::vector<std::set<size_t> >& jacOuterSpar,
                                           const std::vector<std::set<size_t> >& hessOuterSpar,
                                           Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
        ASSERT_TRUE(xNorm.empty() || x.size() == xNorm.size());

        using namespace CppAD;
        using namespace std;
        using CppAD::vector;

        const size_t m = _funOuterNoAtom->Range();

        modelLibOuter->addAtomicFunction(modelLib->asAtomic());

        /**
         * Test zero order
         */
        auto xOrig = makeCGVector(x);

        auto yOrig = _funOuterNoAtom->Forward(0, xOrig);
        auto yOuter = modelLibOuter->ForwardZero(x);

        ASSERT_TRUE(compareValues<double>(yOuter, yOrig, epsilonR, epsilonA));

        /**
         * Jacobian sparsity
         */
        const std::vector<bool> jacSparsityOrig = jacobianSparsityBool<std::vector<bool>, CGD> (*_funOuterNoAtom);

        /**
         * Sparse Jacobian
         */
        CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
        std::vector<size_t> jacOuterRows, jacOuterCols;
        generateSparsityIndexes(jacOuterSpar, jacOuterRows, jacOuterCols);
        std::vector<CGD> jacOrig(jacOuterCols.size());
        _funOuterNoAtom->SparseJacobianReverse(xOrig, jacSparsityOrig,
                                               jacOuterRows, jacOuterCols, jacOrig, work);

        std::vector<double> jacOuter(jacOuterRows.size());
        size_t const* rows, *cols;
        modelLibOuter->SparseJacobian(x, par, jacOuter, &rows, &cols);

        ASSERT_TRUE(compareValues<double>(jacOuter, jacOrig, epsilonR, epsilonA));

        /**
         * Hessian sparsity
         */
        const std::vector<bool> hessSparsityOrig = hessianSparsityBool<std::vector<bool>, CGD> (*_funOuterNoAtom);
        /**
        printSparsityPattern(hessSparsityOrig, "original", n, n);

        std::vector<size_t> hessOuterRows2, hessOuterCols2;
        modelLibOuter->HessianSparsity(hessOuterRows2, hessOuterCols2);
        printSparsityPattern(hessOuterRows2, hessOuterCols2, "outer", n);
        */

        /**
         * Sparse Hessian
         */
        std::vector<CGD> wOrig(m);
        std::vector<double> w(m);
        for (size_t i = 0; i < m; i++) {
            wOrig[i] = 1;
            w[i] = 1;
        }
        CppAD::sparse_hessian_work hessWork;
        std::vector<size_t> hessOuterRows, hessOuterCols;
        generateSparsityIndexes(hessOuterSpar, hessOuterRows, hessOuterCols);
        std::vector<CGD> hessOrig(hessOuterRows.size());
        _funOuterNoAtom->SparseHessian(xOrig, wOrig, hessSparsityOrig,
                                       hessOuterRows, hessOuterCols, hessOrig, hessWork);

        vector<double> hessOuter(hessOuterRows.size());
        modelLibOuter->SparseHessian(x, par, w,
                                     hessOuter, &rows, &cols);

        printVector("hessOuter", hessOuter);
        printVector("hessOrig", hessOrig);

        ASSERT_TRUE(compareValues<double>(hessOuter, hessOrig, epsilonR, epsilonA));
    }

    virtual void prepareAtomicLibAtomicLib(const std::vector<Base>& x,
                                           const std::vector<Base>& par,
                                           const std::vector<Base>& xNorm,
                                           const std::vector<Base>& eqNorm) {
        using CppAD::vector;

        tapeInnerModel(x, par, xNorm, eqNorm);

        const size_t m = _funInner->Range();

        /**
         * Create the dynamic library model
         */
        auto cSourceInner = prepareInnerModelCompilation();

        /**
         * Create the dynamic library
         * (generate and compile source code)
         */
        GccCompiler<double> compiler1;
        prepareTestCompilerFlags(compiler1);
        compiler1.setSourcesFolder("sources_atomiclibatomiclib_" + _modelName);
        compiler1.setSaveToDiskFirst(true);

        ModelLibraryCSourceGen<double> compDynHelp(*cSourceInner);
        compDynHelp.setVerbose(this->verbose_);

        DynamicModelLibraryProcessor<double> p(compDynHelp, "innerModel");
        _dynamicLib = p.createDynamicLibrary(compiler1);
        _modelLib = _dynamicLib->model(_modelName);

        /**
         * Second compiled model
         */
        // independent variables
        std::vector<ADCGD> ax2 = makeADCGVector(x);

        // parameters
        std::vector<ADCGD> ap2 = makeADCGVector(par);

        size_t abort_op_index = 0;
        bool record_compare = true;
        CppAD::Independent(ax2, abort_op_index, record_compare, ap2);

        CGAtomicGenericModel<Base>& innerAtomicFun = _modelLib->asAtomic();
        CGAtomicFun<Base> cgInnerAtomicFun(innerAtomicFun, std::vector<double>(_modelLib->Domain()), _modelLib->Parameters(), true); // required for taping

        std::vector<ADCGD> ay2(m);
        cgInnerAtomicFun(ax2, ay2);
        std::vector<ADCGD> Z2 = modelOuter(ay2, ap2);

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

        DynamicModelLibraryProcessor<double> p2(compDynHelp2, "outerModel");
        GccCompiler<double> compiler2;
        prepareTestCompilerFlags(compiler2);
        compiler2.setSourcesFolder("sources_atomiclibatomiclib_" + _modelName);
        compiler2.setSaveToDiskFirst(true);
        _dynamicLibOuter = p2.createDynamicLibrary(compiler2);

        /**
         * tape the model without atomics
         */
        tapeOuterModelNoAtomics(x, par, xNorm, eqNorm);
    }

    void prepareAtomicLibModelBridge(const std::vector<Base>& x,
                                     const std::vector<Base>& par,
                                     const std::vector<Base>& xNorm,
                                     const std::vector<Base>& eqNorm,
                                     const std::vector<std::set<size_t> >& jacInner = {},
                                     const std::vector<std::set<size_t> >& hessInner = {},
                                     const std::vector<std::set<size_t> >& jacOuter = {},
                                     const std::vector<std::set<size_t> >& hessOuter = {},
                                     bool createOuterReverse2 = true) {

        tapeInnerModel(x, par, xNorm, eqNorm);

        compileInnerAndOuterModelWithAtomicsSameLib(x, par, jacInner, hessInner, jacOuter, hessOuter,
                                                    createOuterReverse2);

        tapeOuterModelNoAtomics(x, par, xNorm, eqNorm);
    }

protected:
    void tapeInnerModel(const std::vector<Base>& x,
                        const std::vector<Base>& par,
                        const std::vector<Base>& xNorm,
                        const std::vector<Base>& eqNorm) {
        const size_t n = x.size();

        // independent variables
        std::vector<ADCGD> ax(n);
        for (size_t j = 0; j < n; j++)
            ax[j] = x[j];

        std::vector<ADCGD> ap;

        size_t abort_op_index = 0;
        bool record_compare = true;
        CppAD::Independent(ax, abort_op_index, record_compare, ap);

        normalizeIndependents(ax, xNorm);

        /**
         * create the CppAD tape as usual
         */
        std::vector<ADCGD> ay = model(ax, ap);

        normalizeDependents(ay, eqNorm);

        // create f: U -> ay and vectors used for derivative calculations
        _funInner = std::make_unique<ADFun<CGD>>();
        _funInner->Dependent(ay);
    }

    void tapeOuterModelNoAtomics(const std::vector<Base>& x,
                                 const std::vector<Base>& par,
                                 const std::vector<Base>& xNorm,
                                 const std::vector<Base>& eqNorm) {
        // independent variables
        std::vector<ADCGD> ax = makeADCGVector(x);

        // parameters
        std::vector<ADCGD> ap = makeADCGVector(par);

        size_t abort_op_index = 0;
        bool record_compare = true;
        CppAD::Independent(ax, abort_op_index, record_compare, ap);

        normalizeIndependents(ax, xNorm);

        /**
         * create the CppAD tape as usual
         */
        std::vector<ADCGD> ay = model(ax, ap);

        normalizeDependents(ay, eqNorm);

        std::vector<ADCGD> az = modelOuter(ay, ap);

        // create f: ax -> az
        _funOuterNoAtom = std::make_unique<ADFun<CGD>>();
        _funOuterNoAtom->Dependent(az);
    }

    template<class Atomic>
    void tapeOuterModelWithAtomic(Atomic atomic,
                                  const std::vector<Base>& x,
                                  const std::vector<Base>& par,
                                  const std::vector<Base>& xNorm,
                                  const std::vector<Base>& eqNorm) {
        const size_t m = _modelLib->Range();

        // independent variables
        std::vector<ADCGD> ax = makeADCGVector(x);

        // parameters
        std::vector<ADCGD> ap = makeADCGVector(par);

        size_t abort_op_index = 0;
        bool record_compare = true;
        CppAD::Independent(ax, abort_op_index, record_compare, ap);

        normalizeIndependents(ax, xNorm);

        /**
         * create the CppAD tape as usual
         */
        std::vector<ADCGD> ay(m);

        // call user function and store CGAtomicLibModel(x) in au[0]
        assert(_modelLib != nullptr);
        assert(_modelLib->getName() ==  _modelName);
        _atomFun = std::make_unique<CGAtomicFun<double>>(_modelLib->asAtomic(), x, _modelLib->Parameters());
        //CGAtomicGenericModel<double>& atomicfun = _modelLib->asAtomic();

        (*_atomFun)(ax, ay); // atomic function

        normalizeDependents(ay, eqNorm);

        std::vector<ADCGD> az = modelOuter(ay, ap);

        // create f: ax -> az
        _funOuterAtom = std::make_unique<ADFun<CGD>>();
        _funOuterAtom->Dependent(az);
    }

    std::unique_ptr<ModelCSourceGen<double>>
    prepareInnerModelCompilation(const std::vector<std::set<size_t> >& jacInner = {},
                                 const std::vector<std::set<size_t> >& hessInner = {}) {
        /**
         * Create the dynamic library model for the inner model
         */
        std::unique_ptr<ModelCSourceGen<double>> cSourceInner(new ModelCSourceGen<double>(*_funInner, _modelName));
        cSourceInner->setCreateForwardZero(true);
        cSourceInner->setCreateForwardOne(true);
        cSourceInner->setCreateReverseOne(true);
        cSourceInner->setCreateReverseTwo(true);
        cSourceInner->setCreateJacobian(true);
        cSourceInner->setCreateHessian(true);
        cSourceInner->setCreateSparseJacobian(true);
        cSourceInner->setCreateSparseHessian(true);

        if (!jacInner.empty()) {
            cSourceInner->setCustomSparseJacobianElements(jacInner);
        }
        if (!hessInner.empty()) {
            cSourceInner->setCustomSparseHessianElements(hessInner);
        }

        return cSourceInner;
    }

    std::unique_ptr<GenericModel<double>> compileInnerModel(const std::string& folderNamePrefix,
                                                            const std::vector<std::set<size_t> >& jacInner = {},
                                                            const std::vector<std::set<size_t> >& hessInner = {}) {

        auto cSourceInner = prepareInnerModelCompilation(jacInner, hessInner);
        /**
         * Create the dynamic library
         * (generate and compile source code)
         */
        ModelLibraryCSourceGen<double> libSrcGen(*cSourceInner);
        libSrcGen.setVerbose(this->verbose_);

        std::string folder = std::string("sources_inner_model") + folderNamePrefix + _modelName;

        DynamicModelLibraryProcessor<double> p(libSrcGen);

        GccCompiler<double> compiler;
        prepareTestCompilerFlags(compiler);
        compiler.setSourcesFolder(folder);
        compiler.setSaveToDiskFirst(true);
        _dynamicLib = p.createDynamicLibrary(compiler);

        return _dynamicLib->model(_modelName);
    }

    void compileInnerAndOuterModelWithAtomicsSameLib(const std::vector<Base>& x,
                                                     const std::vector<Base>& par,
                                                     const std::vector<std::set<size_t> >& jacInner,
                                                     const std::vector<std::set<size_t> >& hessInner,
                                                     const std::vector<std::set<size_t> >& jacOuter,
                                                     const std::vector<std::set<size_t> >& hessOuter,
                                                     bool createOuterReverse2) {
        /*
         * Create the dynamic library model
         */
        auto cSourceInner = prepareInnerModelCompilation(jacInner, hessInner);

        /*
         * Second compiled model
         */
        std::unique_ptr<CGAtomicFunBridge<double>> atomicFun;
        auto cSourceOuter = prepareOuterModelCompilationWithAtomicSameLib(x, par, jacInner, hessInner, jacOuter,
                                                                          hessOuter, createOuterReverse2, atomicFun);

        /*
         * Create the dynamic library
         * (generate and compile source code)
         */
        ModelLibraryCSourceGen<double> libGenInner(*cSourceInner, *cSourceOuter);
        libGenInner.setVerbose(this->verbose_);

        std::string folder = std::string("sources_atomiclibmodelbridge_") + (createOuterReverse2 ? "rev2_" : "dir_") + _modelName;

        DynamicModelLibraryProcessor<double> p(libGenInner);

        GccCompiler<double> compiler;
        prepareTestCompilerFlags(compiler);
        compiler.setSourcesFolder(folder);
        compiler.setSaveToDiskFirst(true);
        _dynamicLib = p.createDynamicLibrary(compiler);
    }

    std::unique_ptr<ModelCSourceGen<double>>
    prepareOuterModelCompilationWithAtomicSameLib(const std::vector<Base>& x,
                                                  const std::vector<Base>& par,
                                                  const std::vector<std::set<size_t> >& jacInner,
                                                  const std::vector<std::set<size_t> >& hessInner,
                                                  const std::vector<std::set<size_t> >& jacOuter,
                                                  const std::vector<std::set<size_t> >& hessOuter,
                                                  bool createOuterReverse2,
                                                  std::unique_ptr<CGAtomicFunBridge<double>>& atomicFun) {
        const size_t m = _funInner->Range();

        // independent variables
        std::vector<ADCGD> axOuter = makeADCGVector(x);

        // parameters
        std::vector<ADCGD> ap2;

        size_t abort_op_index = 0;
        bool record_compare = true;
        CppAD::Independent(axOuter, abort_op_index, record_compare, ap2);

        atomicFun = std::make_unique<CGAtomicFunBridge<double>>(_modelName, *_funInner, true);
        if (!jacInner.empty()) {
            atomicFun->setCustomSparseJacobianElements(jacInner);
        }
        if (!hessInner.empty()) {
            atomicFun->setCustomSparseHessianElements(hessInner);
        }

        std::vector<ADCGD> ayInner(m);
        (*atomicFun)(axOuter, ayInner);
        std::vector<ADCGD> ayOuter = modelOuter(ayInner, ap2);

        // create f: axOuter -> ayOuter
        _funOuterAtom = std::make_unique<ADFun<CGD>>();
        _funOuterAtom->Dependent(ayOuter);

        std::unique_ptr<ModelCSourceGen<double>> cSourceOuter(
                new ModelCSourceGen<double>(*_funOuterAtom, _modelName + "_outer"));

        cSourceOuter->setCreateForwardZero(true);
        cSourceOuter->setCreateForwardOne(true);
        cSourceOuter->setCreateReverseOne(true);
        cSourceOuter->setCreateReverseTwo(createOuterReverse2);
        cSourceOuter->setCreateJacobian(true);
        cSourceOuter->setCreateHessian(true);
        cSourceOuter->setCreateSparseJacobian(true);
        cSourceOuter->setCreateSparseHessian(true);
        if (!jacOuter.empty()) {
            cSourceOuter->setCustomSparseJacobianElements(jacOuter);
        }
        if (!hessOuter.empty()) {
            cSourceOuter->setCustomSparseHessianElements(hessOuter);
        }

        return cSourceOuter;
    }

    static void normalizeIndependents(std::vector<ADCGD>& ax, const std::vector<Base>& xNorm) {
        if (xNorm.empty()) {
            return;
        }

        const size_t n = xNorm.size();
        for (size_t j = 0; j < n; j++) {
            ax[j] *= xNorm[j];
        }
    }

    static void normalizeDependents(std::vector<ADCGD>& ay, const std::vector<Base>& eqNorm) {
        if (eqNorm.empty()) {
            return;
        }
        ASSERT_EQ(ay.size(), eqNorm.size());
        for (size_t i = 0; i < ay.size(); i++)
            ay[i] /= eqNorm[i];

    }
};

} // END cg namespace
} // END CppAD namespace

#endif
