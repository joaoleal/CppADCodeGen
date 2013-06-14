#ifndef CPPADCGOO_TEST_CPPADCGDYNAMICATOMICTEST_INCLUDED
#define CPPADCGOO_TEST_CPPADCGDYNAMICATOMICTEST_INCLUDED
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
        DynamicLib<Base>* _dynamicLib;
    public:

        inline CppADCGDynamicAtomicTest(const std::string& modelName, bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues),
            _modelName(modelName),
            _fun(NULL),
            _dynamicLib(NULL) {
            //this->verbose_ = true;
        }

        virtual std::vector<ADCGD> model(const std::vector<ADCGD>& ind) = 0;

        virtual void TearDown() {
            delete _dynamicLib;
            _dynamicLib = NULL;
            delete _fun;
            _fun = NULL;
        }

        virtual ~CppADCGDynamicAtomicTest() {
            delete _dynamicLib;
            delete _fun;
        }

        void testADFunAtomicLib(const CppAD::vector<Base>& x) {
            CppAD::vector<Base> xNorm(x.size());
            for (size_t i = 0; i < xNorm.size(); i++)
                xNorm[i] = 1.0;
            CppAD::vector<Base> eqNorm;

            testADFunAtomicLib(x, xNorm, eqNorm);
        }

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


            // -----------------------------------------------------------------
            // Free all temporary work space associated with user_atomic objects. 
            // (If there are future calls to user atomic functions, they will 
            // create new temporary work space.)
            CppAD::user_atomic<double>::clear();
        }

    private:

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
            compHelp.setCreateJacobian(false);
            compHelp.setCreateHessian(false);
            compHelp.setCreateSparseJacobian(false);
            compHelp.setCreateSparseHessian(false);

            GccCompiler<double> compiler;

            CLangCompileDynamicHelper<double> compDynHelp(compHelp);
            _dynamicLib = compDynHelp.createDynamicLibrary(compiler);
        }

    };

}

#endif