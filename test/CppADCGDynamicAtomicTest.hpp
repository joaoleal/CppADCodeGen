#ifndef CPPAD_CG_TEST_CPPADCGDYNAMICATOMICTEST_INCLUDED
#define CPPAD_CG_TEST_CPPADCGDYNAMICATOMICTEST_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

#include "CppADCGTest.hpp"

namespace CppAD {

    class CppADCGDynamicAtomicTest : public CppADCGTest {
    public:
        typedef CppADCGTest::Base Base;
        typedef CppADCGTest::CGD CGD;
        typedef CppADCGTest::ADCGD ADCGD;
    protected:
        const std::string _modelName;
        ADFun<CGD>* _fun;
        ADFun<CGD>* _fun2;
        DynamicLib<Base>* _dynamicLib;
        DynamicLib<Base>* _dynamicLib2;
        DynamicLibModel<Base>* _modelLib;
    public:

        inline CppADCGDynamicAtomicTest(const std::string& modelName, bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues),
            _modelName(modelName),
            _fun(NULL),
            _fun2(NULL),
            _dynamicLib(NULL),
            _dynamicLib2(NULL),
            _modelLib(NULL) {
            //this->verbose_ = true;
        }

        virtual std::vector<ADCGD> model(const std::vector<ADCGD>& ind) = 0;

        virtual std::vector<ADCGD> modelOuter(const std::vector<ADCGD>& ind2) {
            std::vector<ADCGD> Z(ind2.size() - 1);

            for (size_t i = 0; i < ind2.size() - 1; i++) {
                Z[i] = 2 * ind2[i];
            }
            Z[Z.size() - 1] += ind2[ind2.size() - 1];

            return Z;
        }

        virtual void TearDown() {
            delete _dynamicLib;
            _dynamicLib = NULL;
            delete _dynamicLib2;
            _dynamicLib2 = NULL;
            delete _modelLib;
            _modelLib = NULL;
            delete _fun;
            _fun = NULL;
            delete _fun2;
            _fun2 = NULL;
        }

        virtual ~CppADCGDynamicAtomicTest() {
            delete _dynamicLib;
            delete _dynamicLib2;
            delete _modelLib;
            delete _fun;
            delete _fun2;
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

            const size_t n = _fun->Domain();
            const size_t m = _fun->Range();

            vector< AD<double> > ax(n);
            for (size_t j = 0; j < n; j++)
                ax[j] = x[j];

            // declare independent variables and start tape recording
            CppAD::Independent(ax);

            vector<AD<double> > ay(m);

            // call user function and store CGAtomicLibModel(x) in au[0] 
            auto_ptr<DynamicLibModel<Base> > modelLib(_dynamicLib->model(_modelName));
            CGAtomicLibModel<double>& atomicfun = modelLib->asAtomic();

            atomicfun(ax, ay);

            // create f2: x -> y and stop tape recording
            ADFun<double> f2(ax, ay);

            /**
             * Test zero order
             */
            vector<CGD> xOrig(n);
            for (size_t j = 0; j < n; j++)
                xOrig[j] = x[j];

            vector<CGD> yOrig = _fun->Forward(0, xOrig);
            vector<double> yInner = modelLib->ForwardZero(x);
            vector<double> yOutter = f2.Forward(0, x);

            compareValues(yInner, yOrig, epsilonR, epsilonA);
            compareValues(yOutter, yOrig, epsilonR, epsilonA);

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

                vector<CGD> y_pOrig = _fun->Forward(1, x_pOrig);
                vector<double> y_pInner = modelLib->ForwardOne(tx);
                vector<double> y_pOutter = f2.Forward(1, x_p);

                x_p[j] = 0;
                x_pOrig[j] = 0;
                tx[j * k1 + 1] = 0;

                compareValues(y_pInner, y_pOrig, epsilonR, epsilonA);
                compareValues(y_pOutter, y_pOrig, epsilonR, epsilonA);
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

                vector<CGD> dwOrig = _fun->Reverse(1, wOrig);
                vector<double> dwInner = modelLib->ReverseOne(tx, ty, w);
                vector<double> dwOutter = f2.Reverse(1, w);

                w[i] = 0;
                wOrig[i] = 0;

                compareValues(dwInner, dwOrig, epsilonR, epsilonA);
                compareValues(dwOutter, dwOrig, epsilonR, epsilonA);
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

                _fun->Forward(1, x_pOrig);
                vector<CGD> dwOrig = _fun->Reverse(2, pyOrig);
                vector<double> dwInner = modelLib->ReverseTwo(tx, ty, py);
                f2.Forward(1, x_p);
                vector<double> dwOutter = f2.Reverse(2, py);

                x_p[j] = 0;
                x_pOrig[j] = 0;
                tx[j * k1 + 1] = 0;

                // only compare second order information
                // (location of the elements is different then if py.size() == m)
                ASSERT_EQ(dwOrig.size(), n * k1);
                ASSERT_EQ(dwOrig.size(), dwInner.size());
                ASSERT_EQ(dwOrig.size(), dwOutter.size());
                for (size_t j = 0; j < n; j++) {
                    nearEqual(dwInner[j * k1], dwOrig[j * k1].getValue());
                    nearEqual(dwOutter[j * k1], dwOrig[j * k1].getValue());
                }
            }

            /**
             * Jacobian
             */
            vector<CGD> jacOrig = _fun->Jacobian(xOrig);
            vector<double> jacOutter = f2.Jacobian(x);
            compareValues(jacOutter, jacOrig, epsilonR, epsilonA);

            /**
             * Jacobian sparsity
             */
            const std::vector<bool> jacSparsityOrig = jacobianForwardSparsity < std::vector<bool>, CGD > (*_fun);
            const std::vector<bool> jacSparsityOutter = jacobianForwardSparsity < std::vector<bool>, double > (f2);

            compareBoolValues(jacSparsityOrig, jacSparsityOutter);

            const std::vector<bool> jacSparsityOrigRev = jacobianReverseSparsity < std::vector<bool>, CGD > (*_fun);
            const std::vector<bool> jacSparsityOutterRev = jacobianReverseSparsity < std::vector<bool>, double > (f2);

            compareBoolValues(jacSparsityOrigRev, jacSparsityOrig);
            compareBoolValues(jacSparsityOrigRev, jacSparsityOutterRev);

            /**
             * Sparse jacobian
             */
            jacOrig = _fun->SparseJacobian(xOrig);
            jacOutter = f2.SparseJacobian(x);
            compareValues(jacOutter, jacOrig, epsilonR, epsilonA);

            // sparse reverse
            std::vector<size_t> row, col;
            generateSparsityIndexes(jacSparsityOutter, m, n, row, col);

            sparse_jacobian_work workOrig;
            jacOrig.resize(row.size());
            _fun->SparseJacobianReverse(xOrig, jacSparsityOrig, row, col, jacOrig, workOrig);

            sparse_jacobian_work work2;
            jacOutter.resize(row.size());
            f2.SparseJacobianReverse(x, jacSparsityOutter, row, col, jacOutter, work2);

            compareValues(jacOutter, jacOrig, epsilonR, epsilonA);

            /**
             * Hessian
             */
            for (size_t i = 0; i < m; i++) {
                w[i] = 1;
                wOrig[i] = 1;
            }
            vector<CGD> hessOrig = _fun->Hessian(xOrig, wOrig);
            vector<double> hessOutter = f2.Hessian(x, w);
            compareValues(hessOutter, hessOrig, epsilonR, epsilonA);

            /**
             * Sparse Hessian
             */
            hessOrig = _fun->SparseHessian(xOrig, wOrig);
            hessOutter = f2.SparseHessian(x, w);
            compareValues(hessOutter, hessOrig, epsilonR, epsilonA);
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
            ASSERT_TRUE(_modelLib != NULL);

            auto_ptr<DynamicLibModel<Base> > modelLibOuter(_dynamicLib2->model(_modelName + "_outer"));
            ASSERT_TRUE(modelLibOuter.get() != NULL);

            test2LevelAtomicLibModel(_modelLib, modelLibOuter.get(),
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

            auto_ptr<DynamicLibModel<Base> > modelLib(_dynamicLib->model(_modelName));
            auto_ptr<DynamicLibModel<Base> > modelLibOuter(_dynamicLib->model(_modelName + "_outer"));

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

            auto_ptr<DynamicLibModel<Base> > modelLib(_dynamicLib->model(_modelName));
            auto_ptr<DynamicLibModel<Base> > modelLibOuter(_dynamicLib->model(_modelName + "_outer"));

            test2LevelAtomicLibModelCustomEls(modelLib.get(), modelLibOuter.get(),
                                              x, xNorm, eqNorm,
                                              jacOuter, hessOuter,
                                              epsilonR, epsilonA);
        }

    private:

        void test2LevelAtomicLibModel(DynamicLibModel<Base>* modelLib,
                                      DynamicLibModel<Base>* modelLibOuter,
                                      const CppAD::vector<Base>& x,
                                      const CppAD::vector<Base>& xNorm,
                                      const CppAD::vector<Base>& eqNorm,
                                      Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
            ASSERT_EQ(x.size(), xNorm.size());

            using namespace CppAD;
            using namespace std;
            using CppAD::vector;

            const size_t n = _fun2->Domain();
            const size_t m = _fun2->Range();

            modelLibOuter->addAtomicFunction(modelLib->asAtomic());


            /**
             * Test zero order
             */
            vector<CGD> xOrig(n);
            for (size_t j = 0; j < n; j++)
                xOrig[j] = x[j];

            vector<CGD> yOrig = _fun2->Forward(0, xOrig);
            vector<double> yOuter = modelLibOuter->ForwardZero(x);

            compareValues(yOuter, yOrig, epsilonR, epsilonA);

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

                vector<CGD> y_pOrig = _fun2->Forward(1, x_pOrig);
                vector<double> y_pOutter = modelLibOuter->ForwardOne(tx);

                x_pOrig[j] = 0;
                tx[j * k1 + 1] = 0;

                compareValues(y_pOutter, y_pOrig, epsilonR, epsilonA);
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

                vector<CGD> dwOrig = _fun2->Reverse(1, wOrig);
                vector<double> dwOutter = modelLibOuter->ReverseOne(tx, ty, w);

                w[i] = 0;
                wOrig[i] = 0;

                compareValues(dwOutter, dwOrig, epsilonR, epsilonA);
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

                _fun2->Forward(1, x_pOrig);
                vector<CGD> dwOrig = _fun2->Reverse(2, pyOrig);
                vector<double> dwOutter = modelLibOuter->ReverseTwo(tx, ty, py);

                x_pOrig[j] = 0;
                tx[j * k1 + 1] = 0;

                // only compare second order information
                // (location of the elements is different then if py.size() == m)
                ASSERT_EQ(dwOrig.size(), n * k1);
                ASSERT_EQ(dwOrig.size(), dwOutter.size());
                for (size_t j = 0; j < n; j++) {
                    nearEqual(dwOutter[j * k1], dwOrig[j * k1].getValue());
                }
            }

            /**
             * Jacobian
             */
            vector<CGD> jacOrig = _fun2->Jacobian(xOrig);
            vector<double> jacOuter = modelLibOuter->Jacobian(x);
            compareValues(jacOuter, jacOrig, epsilonR, epsilonA);

            /**
             * Jacobian sparsity
             */
            const std::vector<bool> jacSparsityOrig = jacobianSparsity < std::vector<bool>, CGD > (*_fun2);
            const std::vector<bool> jacSparsityOuter = modelLibOuter->JacobianSparsityBool();

            compareBoolValues(jacSparsityOrig, jacSparsityOuter);

            /**
             * Sparse jacobian
             */
            jacOrig = _fun2->SparseJacobian(xOrig);
            jacOuter = modelLibOuter->SparseJacobian(x);
            compareValues(jacOuter, jacOrig, epsilonR, epsilonA);

            /**
             * Hessian
             */
            for (size_t i = 0; i < m; i++) {
                w[i] = 1;
                wOrig[i] = 1;
            }
            vector<CGD> hessOrig = _fun2->Hessian(xOrig, wOrig);
            vector<double> hessOuter = modelLibOuter->Hessian(x, w);
            compareValues(hessOuter, hessOrig, epsilonR, epsilonA);

            /**
             * Hessian sparsity
             */
            const std::vector<bool> hessSparsityOrig = hessianSparsity < std::vector<bool>, CGD > (*_fun2);
            const std::vector<bool> hessSparsityOuter = modelLibOuter->HessianSparsityBool();

            compareBoolValues(hessSparsityOrig, hessSparsityOuter);

            /**
             * Sparse Hessian
             */
            hessOrig = _fun2->SparseHessian(xOrig, wOrig);
            hessOuter = modelLibOuter->SparseHessian(x, w);
            compareValues(hessOuter, hessOrig, epsilonR, epsilonA);
        }

        void test2LevelAtomicLibModelCustomEls(DynamicLibModel<Base>* modelLib,
                                               DynamicLibModel<Base>* modelLibOuter,
                                               const CppAD::vector<Base>& x,
                                               const CppAD::vector<Base>& xNorm,
                                               const CppAD::vector<Base>& eqNorm,
                                               const std::vector<std::set<size_t> >& jacOuter,
                                               const std::vector<std::set<size_t> >& hessOuter,
                                               Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
            ASSERT_EQ(x.size(), xNorm.size());

            using namespace CppAD;
            using namespace std;
            using CppAD::vector;

            const size_t n = _fun2->Domain();
            const size_t m = _fun2->Range();

            modelLibOuter->addAtomicFunction(modelLib->asAtomic());

            /**
             * Test zero order
             */
            vector<CGD> xOrig(n);
            for (size_t j = 0; j < n; j++)
                xOrig[j] = x[j];

            vector<CGD> yOrig = _fun2->Forward(0, xOrig);
            vector<double> yOuter = modelLibOuter->ForwardZero(x);

            compareValues(yOuter, yOrig, epsilonR, epsilonA);

            /**
             * Jacobian sparsity
             */
            const std::vector<bool> jacSparsityOrig = jacobianSparsity < std::vector<bool>, CGD > (*_fun2);

            /**
             * Sparse jacobian
             */
            CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
            std::vector<size_t> jacOuterRows, jacOuterCols;
            generateSparsityIndexes(jacOuter, jacOuterRows, jacOuterCols);
            vector<CGD> jacOrig(jacOuterCols.size());
            _fun2->SparseJacobianReverse(xOrig, jacSparsityOrig,
                                         jacOuterRows, jacOuterCols, jacOrig, work);

            std::vector<double> jacOuter(jacOuterRows.size());
            size_t const* rows, *cols;
            modelLibOuter->SparseJacobian(&x[0], x.size(), &jacOuter[0], &rows, &cols, jacOuter.size());

            compareValues(jacOuter, jacOrig, epsilonR, epsilonA);

            /**
             * Hessian sparsity
             */
            const std::vector<bool> hessSparsityOrig = hessianSparsity < std::vector<bool>, CGD > (*_fun2);

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
            generateSparsityIndexes(hessOuter, hessOuterRows, hessOuterCols);
            vector<CGD> hessOrig(hessOuterRows.size());
            _fun2->SparseHessian(xOrig, wOrig, hessSparsityOrig,
                                 hessOuterRows, hessOuterCols, hessOrig, hessWork);

            vector<double> hessOuter(hessOuterRows.size());
            modelLibOuter->SparseHessian(&x[0], x.size(), &w[0], w.size(),
                                         &hessOuter[0], &rows, &cols, hessOuter.size());

            compareValues(hessOuter, hessOrig, epsilonR, epsilonA);
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
            _fun = new ADFun<CGD>();
            _fun->Dependent(Z);

            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            CLangCompileModelHelper<double> compHelp(*_fun, _modelName);

            compHelp.setCreateForwardZero(true);
            compHelp.setCreateForwardOne(true);
            compHelp.setCreateReverseOne(true);
            compHelp.setCreateReverseTwo(true);

            GccCompiler<double> compiler;
            std::vector<std::string> flags;
            flags.push_back("-O0");
            flags.push_back("-g");
            flags.push_back("-ggdb");
            flags.push_back("-D_FORTIFY_SOURCE=2");
            compiler.setCompileFlags(flags);
            compiler.setSourcesFolder("sources_atomiclib_" + _modelName);

            CLangCompileDynamicHelper<double> compDynHelp(compHelp);
            _dynamicLib = compDynHelp.createDynamicLibrary(compiler);
        }

        virtual void prepareAtomicLibAtomicLib(const CppAD::vector<Base>& x,
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
            _fun = new ADFun<CGD>();
            _fun->Dependent(Z);

            /**
             * Create the dynamic library model
             */
            CLangCompileModelHelper<double> compHelp1(*_fun, _modelName);
            compHelp1.setCreateForwardZero(true);
            compHelp1.setCreateForwardOne(true);
            compHelp1.setCreateReverseOne(true);
            compHelp1.setCreateReverseTwo(true);

            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            GccCompiler<double> compiler1;
            std::vector<std::string> flags;
            flags.push_back("-O0");
            flags.push_back("-g");
            flags.push_back("-ggdb");
            flags.push_back("-D_FORTIFY_SOURCE=2");
            compiler1.setCompileFlags(flags);
            compiler1.setSourcesFolder("sources_atomiclibatomiclib_" + _modelName);

            CLangCompileDynamicHelper<double> compDynHelp(compHelp1);
            compDynHelp.setLibraryName("innerModel");
            _dynamicLib = compDynHelp.createDynamicLibrary(compiler1);
            _modelLib = _dynamicLib->model(_modelName);

            /**
             * Second compiled model
             */
            // independent variables
            std::vector<ADCGD> u2(n);
            for (size_t j = 0; j < n; j++)
                u2[j] = x[j];

            CppAD::Independent(u2);

            CGAtomicLibModel<Base>& innerAtomicFun = _modelLib->asAtomic();
            CGAtomicFun<Base> cgInnerAtomicFun(innerAtomicFun, true); // required for taping

            std::vector<ADCGD> ZZ(Z.size());
            cgInnerAtomicFun(u2, ZZ);
            std::vector<ADCGD> Z2 = modelOuter(ZZ);

            // create f: U2 -> Z2
            ADFun<CGD> fun2;
            fun2.Dependent(Z2);

            CLangCompileModelHelper<double> compHelp2(fun2, _modelName + "_outer");
            compHelp2.setCreateForwardZero(true);
            compHelp2.setCreateForwardOne(true);
            compHelp2.setCreateReverseOne(true);
            compHelp2.setCreateReverseTwo(true);
            compHelp2.setCreateJacobian(true);
            compHelp2.setCreateHessian(true);
            compHelp2.setCreateSparseJacobian(true);
            compHelp2.setCreateSparseHessian(true);

            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            GccCompiler<double> compiler2;
            compiler2.setSourcesFolder("sources_atomiclibatomiclib_" + _modelName);
            compiler2.setCompileFlags(flags);

            CLangCompileDynamicHelper<double> compDynHelp2(compHelp2);
            compDynHelp2.setLibraryName("outterModel");
            _dynamicLib2 = compDynHelp2.createDynamicLibrary(compiler2);

            /**
             * 
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
            _fun = new ADFun<CGD>();
            _fun->Dependent(Z);

            /**
             * Create the dynamic library model
             */
            CLangCompileModelHelper<double> compHelp1(*_fun, _modelName);
            if (jacInner.size() > 0) {
                compHelp1.setCustomSparseJacobianElements(jacInner);
            }
            if (hessInner.size() > 0) {
                compHelp1.setCustomSparseHessianElements(hessInner);
            }

            /**
             * Second compiled model
             */
            // independent variables
            std::vector<ADCGD> u2(n);
            for (size_t j = 0; j < n; j++)
                u2[j] = x[j];

            CppAD::Independent(u2);

            CGAtomicFunBridge<double> atomicfun(_modelName, *_fun, true);
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

            CLangCompileModelHelper<double> compHelp2(fun2, _modelName + "_outer");

            compHelp1.setCreateForwardZero(true);
            compHelp1.setCreateForwardOne(true);
            compHelp1.setCreateReverseOne(true);
            compHelp1.setCreateReverseTwo(true);
            //compHelp1.setCreateSparseHessian(true); //not really required

            compHelp2.setCreateForwardZero(true);
            compHelp2.setCreateForwardOne(true);
            compHelp2.setCreateReverseOne(true);
            compHelp2.setCreateReverseTwo(createOuterReverse2);
            compHelp2.setCreateJacobian(true);
            compHelp2.setCreateHessian(true);
            compHelp2.setCreateSparseJacobian(true);
            compHelp2.setCreateSparseHessian(true);
            if (jacOuter.size() > 0) {
                compHelp2.setCustomSparseJacobianElements(jacOuter);
            }
            if (hessOuter.size() > 0) {
                compHelp2.setCustomSparseHessianElements(hessOuter);
            }

            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            GccCompiler<double> compiler;
            std::vector<std::string> flags;
            flags.push_back("-O0");
            flags.push_back("-g");
            flags.push_back("-ggdb");
            flags.push_back("-D_FORTIFY_SOURCE=2");
            compiler.setCompileFlags(flags);
            compiler.setSourcesFolder(std::string("sources_atomiclibmodelbridge_") + (createOuterReverse2 ? "rev2_" : "dir_") + _modelName);

            CLangCompileDynamicHelper<double> compDynHelp(compHelp1);
            compDynHelp.addModel(compHelp2);
            _dynamicLib = compDynHelp.createDynamicLibrary(compiler);

            /**
             * 
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
            std::vector<ADCGD> Z = model(u);
            if (eqNorm.size() > 0) {
                ASSERT_EQ(Z.size(), eqNorm.size());
                for (size_t i = 0; i < Z.size(); i++)
                    Z[i] /= eqNorm[i];
            }

            std::vector<ADCGD> Z2 = modelOuter(Z);

            // create f: U -> Z2
            _fun2 = new ADFun<CGD>();
            _fun2->Dependent(Z2);
        }

    };

}

#endif