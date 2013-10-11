#ifndef CPPAD_CG_TEST_CPPADCGDYNAMICTEST_INCLUDED
#define CPPAD_CG_TEST_CPPADCGDYNAMICTEST_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
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

    class CppADCGDynamicTest : public CppADCGTest {
    public:
        typedef CG<double> CGD;
        typedef AD<CGD> ADCG;
    protected:
        const std::string _name;
    public:

        inline CppADCGDynamicTest(const std::string& testName, bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues),
            _name(testName) {
        }

        virtual std::vector<ADCGD> model(const std::vector<ADCGD>& ind) = 0;

        void testDynamicFull(std::vector<ADCG>& u,
                             const std::vector<double>& x,
                             size_t maxAssignPerFunc = 100,
                             double epsilonR = 1e-14,
                             double epsilonA = 1e-14) {
            const std::vector<double> xNorm(x.size(), 1.0);
            const std::vector<double> eqNorm;

            testDynamicFull(u, x, xNorm, eqNorm, maxAssignPerFunc, epsilonR, epsilonA);
        }

        void testDynamicFull(std::vector<ADCG>& u,
                             const std::vector<double>& x,
                             const std::vector<double>& xNorm,
                             const std::vector<double>& eqNorm,
                             size_t maxAssignPerFunc = 100,
                             double epsilonR = 1e-14,
                             double epsilonA = 1e-14) {
            ASSERT_EQ(u.size(), x.size());
            ASSERT_EQ(x.size(), xNorm.size());

            using namespace std;

            // use a special object for source code generation
            CppAD::Independent(u);

            for (size_t i = 0; i < u.size(); i++)
                u[i] *= xNorm[i];

            // dependent variable vector 
            std::vector<ADCG> Z = model(u);

            if (eqNorm.size() > 0) {
                ASSERT_EQ(Z.size(), eqNorm.size());
                for (size_t i = 0; i < Z.size(); i++)
                    Z[i] /= eqNorm[i];
            }

            /**
             * create the CppAD tape as usual
             */
            // create f: U -> Z and vectors used for derivative calculations
            ADFun<CGD> fun;
            fun.Dependent(Z);

            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            CLangCompileModelHelper<double> compHelp(fun, _name + "dynamic");

            compHelp.setCreateForwardZero(true);
            compHelp.setCreateJacobian(true);
            compHelp.setCreateHessian(true);
            compHelp.setCreateSparseJacobian(true);
            compHelp.setCreateSparseHessian(true);
            compHelp.setCreateForwardOne(true);
            compHelp.setCreateReverseOne(true);
            compHelp.setCreateReverseTwo(true);
            compHelp.setMaxAssignmentsPerFunc(maxAssignPerFunc);

            GccCompiler<double> compiler;
            compiler.setSourcesFolder("sources_" + _name + "_1");

            CLangCompileDynamicHelper<double> compDynHelp(compHelp);
            DynamicLib<double>* dynamicLib = compDynHelp.createDynamicLibrary(compiler);

            /**
             * test the library
             */
            DynamicLibModel<double>* model = dynamicLib->model(_name + "dynamic");
            ASSERT_TRUE(model != NULL);

            // dimensions
            ASSERT_EQ(model->Domain(), fun.Domain());
            ASSERT_EQ(model->Range(), fun.Range());

            /**
             */
            std::vector<CGD> x2(x.size());
            for (size_t i = 0; i < x.size(); i++) {
                x2[i] = x[i];
            }

            // forward zero
            std::vector<CGD> dep = fun.Forward(0, x2);
            std::vector<double> depCGen = model->ForwardZero(x);
            ASSERT_TRUE(compareValues(depCGen, dep, epsilonR, epsilonA));

            // Jacobian
            std::vector<CGD> jac = fun.Jacobian(x2);
            depCGen = model->Jacobian(x);
            ASSERT_TRUE(compareValues(depCGen, jac, epsilonR, epsilonA));

            // Hessian
            std::vector<CGD> w2(Z.size(), 1.0);
            std::vector<double> w(Z.size(), 1.0);

            std::vector<CGD> hess = fun.Hessian(x2, w2);
            depCGen = model->Hessian(x, w);
            ASSERT_TRUE(compareValues(depCGen, hess, epsilonR, epsilonA));

            // sparse Jacobian
            std::vector<double> jacCGen;
            std::vector<size_t> row, col;
            model->SparseJacobian(x, jacCGen, row, col);
            std::vector<double> jacCGenDense(jac.size());
            for (size_t i = 0; i < jacCGen.size(); i++) {
                jacCGenDense[row[i] * x.size() + col[i]] = jacCGen[i];
            }

            ASSERT_TRUE(compareValues(jacCGenDense, jac, epsilonR, epsilonA));

            // sparse Hessian
            std::vector<double> hessCGen;
            model->SparseHessian(x, w, hessCGen, row, col);
            std::vector<double> hessCGenDense(hess.size());
            for (size_t i = 0; i < hessCGen.size(); i++) {
                hessCGenDense[row[i] * x.size() + col[i]] = hessCGen[i];
            }

            ASSERT_TRUE(compareValues(hessCGenDense, hess, epsilonR, epsilonA));

            delete model;
            delete dynamicLib;
        }

        void testDynamicCustomElements(std::vector<ADCG>& u,
                                       const std::vector<double>& x,
                                       const std::vector<size_t>& jacRow, const std::vector<size_t>& jacCol,
                                       const std::vector<size_t>& hessRow, const std::vector<size_t>& hessCol) {
            using namespace std;

            CppAD::Independent(u);

            // dependent variable vector 
            std::vector<ADCG> Z = model(u);

            /**
             * create the CppAD tape as usual
             */
            // create f: U -> Z and vectors used for derivative calculations
            ADFun<CGD> fun(u, Z);

            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            CLangCompileModelHelper<double> compHelp(fun, _name + "dynamic2");

            compHelp.setCreateSparseJacobian(true);
            compHelp.setCustomSparseJacobianElements(jacRow, jacCol);

            compHelp.setCreateSparseHessian(true);
            compHelp.setCustomSparseHessianElements(hessRow, hessCol);

            GccCompiler<double> compiler;
            compiler.setSourcesFolder("sources_" + _name + "_2");

            CLangCompileDynamicHelper<double> compDynHelp(compHelp);
            compDynHelp.setLibraryName("cppad_cg_model_2");
            DynamicLib<double>* dynamicLib = compDynHelp.createDynamicLibrary(compiler);

            /**
             * test the library
             */
            DynamicLibModel<double>* model = dynamicLib->model(_name + "dynamic2");
            ASSERT_TRUE(model != NULL);

            // dimensions
            ASSERT_EQ(model->Domain(), fun.Domain());
            ASSERT_EQ(model->Range(), fun.Range());

            /**
             */
            std::vector<CGD> x2(x.size());
            for (size_t i = 0; i < x.size(); i++) {
                x2[i] = x[i];
            }

            std::vector<size_t> row, col;

            // sparse Jacobian
            std::vector<double> jacCGen;
            model->SparseJacobian(x, jacCGen, row, col);
            std::vector<CG<double> > jacSparse(row.size());

            std::vector<CGD> jac = fun.Jacobian(x2);
            for (size_t i = 0; i < row.size(); i++) {
                jacSparse[i] = jac[row[i] * x.size() + col[i]];
            }

            ASSERT_TRUE(compareValues(jacCGen, jacSparse));

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

            ASSERT_TRUE(compareValues(hessCGen, hessSparse));

            delete model;
            delete dynamicLib;
        }

        inline ::testing::AssertionResult compareValues(const std::vector<double>& depCGen,
                                                        const std::vector<CppAD::CG<double> >& dep,
                                                        double epsilonR = 1e-14, double epsilonA = 1e-14) {

            std::vector<double> depd(dep.size());

            for (size_t i = 0; i < depd.size(); i++) {
                depd[i] = dep[i].getValue();
            }

            return CppADCGTest::compareValues(depCGen, depd, epsilonR, epsilonA);
        }

    };

}
#endif