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

namespace CppAD {

    class CppADCGDynamicAtomicNestedTest : public CppADCGTest {
    public:
        typedef CppADCGTest::Base Base;
        typedef CppADCGTest::CGD CGD;
        typedef CppADCGTest::ADCGD ADCGD;
    protected:
        const std::string _modelName;
        ADFun<CGD>* _fun;
        ADFun<CGD>* _fun2;
        checkpoint<CGD>* _atomicInnerModel;
        DynamicLib<Base>* _dynamicLib;
        DynamicLib<Base>* _dynamicLib2;
        GenericModel<Base>* _modelLib;
    public:

        inline CppADCGDynamicAtomicNestedTest(const std::string& modelName,
                                              bool verbose = false,
                                              bool printValues = false) :
            CppADCGTest(verbose, printValues),
            _modelName(modelName),
            _fun(nullptr),
            _fun2(nullptr),
            _atomicInnerModel(nullptr),
            _dynamicLib(nullptr),
            _dynamicLib2(nullptr),
            _modelLib(nullptr) {
            //this->verbose_ = true;
        }

        virtual std::vector<ADCGD> modelInner(const std::vector<ADCGD>& xInner) = 0;

        virtual std::vector<ADCGD> modelOuter(const std::vector<ADCGD>& xOuter,
                                              atomic_base<CGD>& atomicInnerModel,
                                              size_t xInnerSize,
                                              size_t yInnerSize) {
            std::vector<ADCGD> yOuter(yInnerSize * 2);
            std::vector<ADCGD> xInner(xInnerSize), yInner(yInnerSize);

            std::copy(xOuter.begin(), xOuter.begin() + xInnerSize, xInner.begin());
            atomicInnerModel(xInner, yInner);
            std::copy(yInner.begin(), yInner.end(), yOuter.begin());

            std::copy(xOuter.begin() + xInnerSize, xOuter.end(), xInner.begin());
            atomicInnerModel(xInner, yInner);
            std::copy(yInner.begin(), yInner.end(), yOuter.begin() + yInner.size());

            return yOuter;
        }

        virtual void TearDown() {
            delete _dynamicLib;
            _dynamicLib = nullptr;
            delete _dynamicLib2;
            _dynamicLib2 = nullptr;
            delete _modelLib;
            _modelLib = nullptr;
            delete _atomicInnerModel;
            _atomicInnerModel = nullptr;
            delete _fun;
            _fun = nullptr;
            delete _fun2;
            _fun2 = nullptr;
        }

        virtual ~CppADCGDynamicAtomicNestedTest() {
            delete _dynamicLib;
            delete _dynamicLib2;
            delete _modelLib;
            delete _atomicInnerModel;
            delete _fun;
            delete _fun2;
        }

        /**
         * Test 2 models in 2 dynamic libraries
         */
        void testAtomicLibAtomicLib(const CppAD::vector<Base>& xOuter,
                                    const CppAD::vector<Base>& xInner,
                                    const CppAD::vector<Base>& xNorm,
                                    const CppAD::vector<Base>& eqNorm,
                                    Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
            using namespace std;

            prepareAtomicLibAtomicLib(xOuter, xInner, xNorm, eqNorm);
            ASSERT_TRUE(_modelLib != nullptr);

            unique_ptr<GenericModel<Base> > modelLibOuter(_dynamicLib2->model(_modelName + "_outer"));
            ASSERT_TRUE(modelLibOuter.get() != nullptr);

            test2LevelAtomicLibModel(_modelLib, modelLibOuter.get(),
                                     xOuter, xInner, epsilonR, epsilonA);
        }

        /**
         * Test 2 models in the same dynamic library
         */
        void testAtomicLibModelBridge(const CppAD::vector<Base>& xOuter,
                                      const CppAD::vector<Base>& xInner,
                                      const CppAD::vector<Base>& xNorm,
                                      const CppAD::vector<Base>& eqNorm,
                                      Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
            using namespace std;

            prepareAtomicLibModelBridge(xOuter, xInner, xNorm, eqNorm);

            unique_ptr<GenericModel<Base> > modelLib(_dynamicLib->model(_modelName));
            unique_ptr<GenericModel<Base> > modelLibOuter(_dynamicLib->model(_modelName + "_outer"));

            test2LevelAtomicLibModel(modelLib.get(), modelLibOuter.get(),
                                     xOuter,
                                     xInner, epsilonR, epsilonA);
        }

        /**
         * Test 2 models in the same dynamic library
         */
        void testAtomicLibModelBridgeCustom(const CppAD::vector<Base>& xOuter,
                                            const CppAD::vector<Base>& xInner,
                                            const CppAD::vector<Base>& xNorm,
                                            const CppAD::vector<Base>& eqNorm,
                                            const std::vector<std::set<size_t> >& jacInner,
                                            const std::vector<std::set<size_t> >& hessInner,
                                            const std::vector<std::set<size_t> >& jacOuter,
                                            const std::vector<std::set<size_t> >& hessOuter,
                                            bool createOuterReverse2,
                                            Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
            using namespace std;

            prepareAtomicLibModelBridge(xOuter,
                                        xInner, xNorm, eqNorm,
                                        jacInner, hessInner,
                                        jacOuter, hessOuter,
                                        createOuterReverse2);

            unique_ptr<GenericModel<Base> > modelLib(_dynamicLib->model(_modelName));
            unique_ptr<GenericModel<Base> > modelLibOuter(_dynamicLib->model(_modelName + "_outer"));

            test2LevelAtomicLibModelCustomEls(modelLib.get(), modelLibOuter.get(),
                                              xOuter,
                                              jacOuter, hessOuter,
                                              epsilonR, epsilonA);
        }

    private:

        void test2LevelAtomicLibModel(GenericModel<Base>* modelLib,
                                      GenericModel<Base>* modelLibOuter,
                                      const CppAD::vector<Base>& xOuter,
                                      const CppAD::vector<Base>& xInner,
                                      Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
            using namespace CppAD;
            using namespace CppAD::extra;
            using namespace std;
            using CppAD::vector;

            const size_t n = _fun2->Domain();
            const size_t m = _fun2->Range();

            modelLibOuter->addExternalModel(*modelLib);


            /**
             * Test zero order
             */
            vector<CGD> xOrig(n);
            for (size_t j = 0; j < n; j++)
                xOrig[j] = xOuter[j];

            vector<CGD> yOrig = _fun2->Forward(0, xOrig);
            vector<double> yOuter = modelLibOuter->ForwardZero(xOuter);

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
                tx[j * k1] = xOuter[j]; // zero order
            for (size_t j = 0; j < n; j++)
                tx[j * k1 + 1] = 0; // first order

            for (size_t j = 0; j < n; j++) {
                x_pOrig[j] = 1;
                tx[j * k1 + 1] = 1;

                vector<CGD> y_pOrig = _fun2->Forward(1, x_pOrig);
                vector<double> y_pOutter = modelLibOuter->ForwardOne(tx);

                x_pOrig[j] = 0;
                tx[j * k1 + 1] = 0;

                ASSERT_TRUE(compareValues(y_pOutter, y_pOrig, epsilonR, epsilonA));
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
                tx[j * k1] = xOuter[j]; // zero order
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

                ASSERT_TRUE(compareValues(dwOutter, dwOrig, epsilonR, epsilonA));
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
                tx[j * k1] = xOuter[j]; // zero order
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
                    ASSERT_TRUE(nearEqual(dwOutter[j * k1], dwOrig[j * k1].getValue()));
                }
            }

            /**
             * Jacobian
             */
            vector<CGD> jacOrig = _fun2->Jacobian(xOrig);
            vector<double> jacOuter = modelLibOuter->Jacobian(xOuter);
            ASSERT_TRUE(compareValues(jacOuter, jacOrig, epsilonR, epsilonA));

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
            jacOuter = modelLibOuter->SparseJacobian(xOuter);
            ASSERT_TRUE(compareValues(jacOuter, jacOrig, epsilonR, epsilonA));

            /**
             * Hessian
             */
            for (size_t i = 0; i < m; i++) {
                w[i] = 1;
                wOrig[i] = 1;
            }
            vector<CGD> hessOrig = _fun2->Hessian(xOrig, wOrig);
            vector<double> hessOuter = modelLibOuter->Hessian(xOuter, w);
            ASSERT_TRUE(compareValues(hessOuter, hessOrig, epsilonR, epsilonA));

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
            hessOuter = modelLibOuter->SparseHessian(xOuter, w);
            ASSERT_TRUE(compareValues(hessOuter, hessOrig, epsilonR, epsilonA));
        }

        void test2LevelAtomicLibModelCustomEls(GenericModel<Base>* modelLib,
                                               GenericModel<Base>* modelLibOuter,
                                               const CppAD::vector<Base>& xOuter,
                                               const std::vector<std::set<size_t> >& jacOuterSparsity,
                                               const std::vector<std::set<size_t> >& hessOuterSparsity,
                                               Base epsilonR = 1e-14, Base epsilonA = 1e-14) {
            using namespace CppAD;
            using namespace CppAD::extra;
            using namespace std;
            using CppAD::vector;

            const size_t n2 = _fun2->Domain();
            const size_t m2 = _fun2->Range();

            modelLibOuter->addExternalModel(*modelLib);

            /**
             * Test zero order
             */
            vector<CGD> xOrig(n2);
            for (size_t j = 0; j < n2; j++)
                xOrig[j] = xOuter[j];

            vector<CGD> yOrig = _fun2->Forward(0, xOrig);
            vector<double> yOuter = modelLibOuter->ForwardZero(xOuter);

            ASSERT_TRUE(compareValues(yOuter, yOrig, epsilonR, epsilonA));

            /**
             * Jacobian sparsity
             */
            const std::vector<bool> jacSparsityOrig = jacobianSparsity < std::vector<bool>, CGD > (*_fun2);

            /**
             * Sparse jacobian
             */
            CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
            std::vector<size_t> jacOuterRows, jacOuterCols;
            generateSparsityIndexes(jacOuterSparsity, jacOuterRows, jacOuterCols);
            vector<CGD> jacOrig(jacOuterCols.size());
            _fun2->SparseJacobianReverse(xOrig, jacSparsityOrig,
                                         jacOuterRows, jacOuterCols, jacOrig, work);

            std::vector<double> jacOuter(jacOuterRows.size());
            size_t const* rows, *cols;
            modelLibOuter->SparseJacobian(&xOuter[0], xOuter.size(), &jacOuter[0], &rows, &cols, jacOuter.size());

            ASSERT_TRUE(compareValues(jacOuter, jacOrig, epsilonR, epsilonA));

            /**
             * Hessian sparsity
             */
            const std::vector<bool> hessSparsityOrig = hessianSparsity < std::vector<bool>, CGD > (*_fun2);

            /**
             * Sparse Hessian
             */
            vector<CGD> wOrig(m2);
            vector<double> w(m2);
            for (size_t i = 0; i < m2; i++) {
                wOrig[i] = 1;
                w[i] = 1;
            }
            CppAD::sparse_hessian_work hessWork;
            std::vector<size_t> hessOuterRows, hessOuterCols;
            generateSparsityIndexes(hessOuterSparsity, hessOuterRows, hessOuterCols);
            vector<CGD> hessOrig(hessOuterRows.size());
            _fun2->SparseHessian(xOrig, wOrig, hessSparsityOrig,
                                 hessOuterRows, hessOuterCols, hessOrig, hessWork);

            vector<double> hessOuter(hessOuterRows.size());
            modelLibOuter->SparseHessian(&xOuter[0], xOuter.size(), &w[0], w.size(),
                                         &hessOuter[0], &rows, &cols, hessOuter.size());

            ASSERT_TRUE(compareValues(hessOuter, hessOrig, epsilonR, epsilonA));
        }

        virtual void prepareAtomicLibAtomicLib(const CppAD::vector<Base>& xOuter,
                                               const CppAD::vector<Base>& xInner,
                                               const CppAD::vector<Base>& xInnerNorm,
                                               const CppAD::vector<Base>& eqInnerNorm) {
            const size_t n = xInner.size();
            const size_t m = eqInnerNorm.size();

            // independent variables
            std::vector<ADCGD> u(n);
            for (size_t j = 0; j < n; j++)
                u[j] = xInner[j];

            CppAD::Independent(u);

            /**
             * create the CppAD tape as usual
             */
            struct InnerModelStruct innerModel(*this, xInnerNorm, eqInnerNorm);
            std::vector<ADCGD> Z(m);
            innerModel(u, Z);

            // create f: U -> Z and vectors used for derivative calculations
            _fun = new ADFun<CGD>();
            _fun->Dependent(Z);

            /**
             * Create the dynamic library model
             */
            ModelCSourceGen<double> compHelp1(*_fun, _modelName);
            compHelp1.setCreateForwardZero(true);
            compHelp1.setCreateForwardOne(true);
            compHelp1.setCreateReverseOne(true);
            compHelp1.setCreateReverseTwo(true);

            /**
             * generate source code
             */
            ModelLibraryCSourceGen<double> compDynHelp(compHelp1);

            SaveFilesModelLibraryProcessor<double>::saveLibrarySourcesTo(compDynHelp, "nested_sources_atomiclibatomiclib_" + _modelName);

            /**
             * Create the dynamic library
             * (compile source code)
             */
            DynamicModelLibraryProcessor<double> p(compDynHelp, "innerModel");

            GccCompiler<double> compiler1;
            std::vector<std::string> flags;
            flags.push_back("-O0");
            flags.push_back("-g");
            flags.push_back("-ggdb");
            flags.push_back("-D_FORTIFY_SOURCE=2");
            compiler1.setCompileFlags(flags);

            _dynamicLib = p.createDynamicLibrary(compiler1);
            _modelLib = _dynamicLib->model(_modelName);

            /**
             * Second compiled model
             */
            // independent variables
            const size_t n2 = xOuter.size();
            std::vector<ADCGD> u2(n2);
            for (size_t j = 0; j < n2; j++)
                u2[j] = xOuter[j];

            CppAD::Independent(u2);

            CGAtomicGenericModel<Base>& innerAtomicFun = _modelLib->asAtomic();
            CGAtomicFun<Base> cgInnerAtomicFun(innerAtomicFun, true); // required for taping

            std::vector<ADCGD> Z2 = modelOuter(u2, cgInnerAtomicFun, n, m);

            // create f: U2 -> Z2
            ADFun<CGD> fun2;
            fun2.Dependent(Z2);

            ModelCSourceGen<double> compHelp2(fun2, _modelName + "_outer");
            compHelp2.setCreateForwardZero(true);
            compHelp2.setCreateForwardOne(true);
            compHelp2.setCreateReverseOne(true);
            compHelp2.setCreateReverseTwo(true);
            compHelp2.setCreateJacobian(true);
            compHelp2.setCreateHessian(true);
            compHelp2.setCreateSparseJacobian(true);
            compHelp2.setCreateSparseHessian(true);

            /**
             * generate source code
             */
            ModelLibraryCSourceGen<double> compDynHelp2(compHelp2);
            SaveFilesModelLibraryProcessor<double>::saveLibrarySourcesTo(compDynHelp, "nested_sources_atomiclibatomiclib_" + _modelName);

            /**
             * Create the dynamic library
             * (compile source code)
             */
            GccCompiler<double> compiler2;
            compiler2.setCompileFlags(flags);

            DynamicModelLibraryProcessor<double> p2(compDynHelp2, "outterModel");
            _dynamicLib2 = p2.createDynamicLibrary(compiler2);

            /**
             * 
             */
            tapeOuterModel(xOuter, xInner, xInnerNorm, eqInnerNorm);
        }

        virtual void prepareAtomicLibModelBridge(const CppAD::vector<Base>& xOuter,
                                                 const CppAD::vector<Base>& xInner,
                                                 const CppAD::vector<Base>& xNorm,
                                                 const CppAD::vector<Base>& eqNorm) {
            std::vector<std::set<size_t> > jacInner, hessInner;
            std::vector<std::set<size_t> > jacOuter, hessOuter;
            prepareAtomicLibModelBridge(xOuter, xInner, xNorm, eqNorm,
                                        jacInner, hessInner,
                                        jacOuter, hessOuter,
                                        true);
        }

        virtual void prepareAtomicLibModelBridge(const CppAD::vector<Base>& xOuter,
                                                 const CppAD::vector<Base>& xInner,
                                                 const CppAD::vector<Base>& xInnerNorm,
                                                 const CppAD::vector<Base>& eqInnerNorm,
                                                 const std::vector<std::set<size_t> >& jacInner,
                                                 const std::vector<std::set<size_t> >& hessInner,
                                                 const std::vector<std::set<size_t> >& jacOuter,
                                                 const std::vector<std::set<size_t> >& hessOuter,
                                                 bool createOuterReverse2) {
            const size_t n = xInner.size();
            const size_t m = eqInnerNorm.size();

            // independent variables
            std::vector<ADCGD> u(n);
            for (size_t j = 0; j < n; j++)
                u[j] = xInner[j];

            CppAD::Independent(u);

            /**
             * create the CppAD tape as usual
             */
            struct InnerModelStruct innerModel(*this, xInnerNorm, eqInnerNorm);
            std::vector<ADCGD> Z(eqInnerNorm.size());
            innerModel(u, Z);

            // create f: U -> Z and vectors used for derivative calculations
            _fun = new ADFun<CGD>();
            _fun->Dependent(Z);

            /**
             * Create the dynamic library model
             */
            ModelCSourceGen<double> compHelp1(*_fun, _modelName);
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
            std::vector<ADCGD> u2(xOuter.size());
            for (size_t j = 0; j < u2.size(); j++)
                u2[j] = xOuter[j];

            CppAD::Independent(u2);

            CGAtomicFunBridge<double> atomicfun(_modelName, *_fun, true);
            if (jacInner.size() > 0) {
                atomicfun.setCustomSparseJacobianElements(jacInner);
            }
            if (hessInner.size() > 0) {
                atomicfun.setCustomSparseHessianElements(hessInner);
            }

            std::vector<ADCGD> Z2 = modelOuter(u2, atomicfun, n, m);

            // create f: U2 -> Z2
            ADFun<CGD> fun2;
            fun2.Dependent(Z2);

            ModelCSourceGen<double> compHelp2(fun2, _modelName + "_outer");

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
             * generate source code
             */
            ModelLibraryCSourceGen<double> compDynHelp(compHelp1);
            compDynHelp.addModel(compHelp2);
            std::string folder = std::string("nested_sources_atomiclibmodelbridge_") + (createOuterReverse2 ? "rev2_" : "dir_") + _modelName;
            SaveFilesModelLibraryProcessor<double>::saveLibrarySourcesTo(compDynHelp, folder);

            /**
             * Create the dynamic library
             * (compile source code)
             */
            GccCompiler<double> compiler;
            std::vector<std::string> flags;
            flags.push_back("-O0");
            flags.push_back("-g");
            flags.push_back("-ggdb");
            flags.push_back("-D_FORTIFY_SOURCE=2");
            compiler.setCompileFlags(flags);

            DynamicModelLibraryProcessor<double> p(compDynHelp);

            _dynamicLib = p.createDynamicLibrary(compiler);

            /**
             * 
             */
            tapeOuterModel(xOuter, xInner, xInnerNorm, eqInnerNorm);
        }

        struct InnerModelStruct {
            CppADCGDynamicAtomicNestedTest& tester;
            const CppAD::vector<Base>& xNorm;
            const CppAD::vector<Base>& eqNorm;

            InnerModelStruct(CppADCGDynamicAtomicNestedTest& t,
                             const CppAD::vector<Base>& xN,
                             const CppAD::vector<Base>& eqN) :
                tester(t),
                xNorm(xN),
                eqNorm(eqN) {
            }

            void operator()(const std::vector<ADCGD>& ax, std::vector<ADCGD>& ay) {
                std::vector<ADCGD> u(ax.size());
                for (size_t j = 0; j < ax.size(); j++)
                    u[j] = ax[j] * xNorm[j];

                ay = tester.modelInner(ax);

                if (eqNorm.size() > 0) {
                    ASSERT_EQ(ay.size(), eqNorm.size());
                    for (size_t i = 0; i < ay.size(); i++)
                        ay[i] /= eqNorm[i];
                }
            }
        };

        void tapeOuterModel(const CppAD::vector<Base>& xOuter,
                            const CppAD::vector<Base>& xInner,
                            const CppAD::vector<Base>& xInnerNorm,
                            const CppAD::vector<Base>& eqInnerNorm) {
            const size_t n = xInner.size();
            const size_t m = eqInnerNorm.size();

            // independent variables
            std::vector<ADCGD> u(n);
            for (size_t j = 0; j < n; j++)
                u[j] = xInner[j];

            // dependent
            std::vector<ADCGD> Z(eqInnerNorm.size());

            //atomic inner model
            struct InnerModelStruct innerModel(*this, xInnerNorm, eqInnerNorm);
            _atomicInnerModel = new checkpoint<CGD>("", innerModel, u, Z);

            /**
             * create the CppAD tape as usual
             */
            size_t n2 = xOuter.size();
            std::vector<ADCGD> u2(n2);
            for (size_t j = 0; j < n2; j++)
                u2[j] = xOuter[j];

            CppAD::Independent(u2);

            std::vector<ADCGD> Z2 = modelOuter(u2, *_atomicInnerModel, n, m);

            // create f: U -> Z2
            _fun2 = new ADFun<CGD>();
            _fun2->Dependent(Z2);
        }

    };

}

#endif