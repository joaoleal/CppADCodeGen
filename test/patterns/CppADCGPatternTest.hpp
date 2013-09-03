#ifndef CPPAD_CG_TEST_CPPADCGPATTERNTEST_INCLUDED
#define	CPPAD_CG_TEST_CPPADCGPATTERNTEST_INCLUDED
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

    class CppADCGPatternTest : public CppADCGTest {
    public:
        typedef double Base;
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCGD;
    private:

        enum TEST_TYPE {
            MUST_PASS, MUST_FAIL, IGNORE
        };
    public:

        inline CppADCGPatternTest(bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues) {
        }

        void testPatternDetection(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat),
                                  size_t m,
                                  size_t n,
                                  size_t repeat,
                                  size_t n_loops = 1) {
            using namespace CppAD;

            //size_t m2 = repeat * m;
            size_t n2 = repeat * n;

            /**
             * Tape model
             */
            std::vector<ADCGD> x(n2);
            for (size_t j = 0; j < n2; j++)
                x[j] = 0.5;
            CppAD::Independent(x);

            std::vector<ADCGD> y = (*model)(x, repeat);

            ADFun<CGD> fun;
            fun.Dependent(y);

            testResults(fun, m, repeat, n_loops);
        }

        void testLibCreation(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat),
                             size_t m,
                             size_t n,
                             size_t repeat,
                             const std::string& libName,
                             bool jacobian = true,
                             bool hessian = true) {
            using namespace CppAD;

            //size_t m2 = repeat * m;
            size_t n2 = repeat * n;

            /**
             * Tape model
             */
            std::vector<ADCGD> x(n2);
            for (size_t j = 0; j < n2; j++)
                x[j] = 0.5;
            CppAD::Independent(x);

            std::vector<ADCGD> y = (*model)(x, repeat);

            ADFun<CGD> fun;
            fun.Dependent(y);

            testSourceCodeGen(fun, m, repeat, libName, FORWARD, jacobian ? MUST_PASS : MUST_FAIL, hessian ? MUST_PASS : MUST_FAIL);
            if (jacobian) {
                testSourceCodeGen(fun, m, repeat, libName, REVERSE, MUST_PASS, IGNORE);
            }
            if (hessian) {
                testSourceCodeGen(fun, m, repeat, libName, FORWARD, IGNORE, MUST_PASS, true);
            }
        }

        void testPatternDetectionWithAtomics(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<Base>*>& atoms),
                                             const std::vector<atomic_base<Base>* >& atoms,
                                             size_t m,
                                             size_t n,
                                             size_t repeat,
                                             size_t n_loops = 1) {
            using namespace CppAD;

            std::vector<CGAbstractAtomicFun<double>*> atomics(atoms.size());
            for (size_t a = 0; a < atoms.size(); a++) {
                atomics[a] = new CGAtomicFun<Base>(*atoms[a], true);
            }

            //size_t m2 = repeat * m;
            size_t n2 = repeat * n;

            /**
             * Tape model
             */
            std::vector<ADCGD> x(n2);
            for (size_t j = 0; j < n2; j++)
                x[j] = 0.5;
            CppAD::Independent(x);

            std::vector<ADCGD> y = (*model)(x, repeat, atomics);

            ADFun<CGD> fun;
            fun.Dependent(y);

            testResults(fun, m, repeat, n_loops);

            for (size_t a = 0; a < atomics.size(); a++) {
                delete atomics[a];
            }
        }

        void testLibCreationWithAtomics(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<Base>*>& atoms),
                                        const std::vector<atomic_base<Base>* >& atoms,
                                        size_t m,
                                        size_t n,
                                        size_t repeat,
                                        const std::string& name) {
            using namespace CppAD;

            std::vector<CGAbstractAtomicFun<double>*> atomics(atoms.size());
            for (size_t a = 0; a < atoms.size(); a++) {
                atomics[a] = new CGAtomicFun<Base>(*atoms[a], true);
            }

            //size_t m2 = repeat * m;
            size_t n2 = repeat * n;

            /**
             * Tape model
             */
            std::vector<ADCGD> x(n2);
            for (size_t j = 0; j < n2; j++)
                x[j] = 0.5;
            CppAD::Independent(x);

            std::vector<ADCGD> y = (*model)(x, repeat, atomics);

            ADFun<CGD> fun;
            fun.Dependent(y);

            testSourceCodeGen(fun, m, repeat, name, atoms, FORWARD);
            testSourceCodeGen(fun, m, repeat, name, atoms, REVERSE);
            testSourceCodeGen(fun, m, repeat, name, atoms, FORWARD, IGNORE, MUST_PASS, true);

            for (size_t a = 0; a < atomics.size(); a++) {
                delete atomics[a];
            }
        }

    private:

        std::vector<std::set<size_t> > createRelatedDepCandidates(size_t m, size_t repeat) {
            std::vector<std::set<size_t> > relatedDepCandidates(m);
            for (size_t i = 0; i < repeat; i++) {
                for (size_t ii = 0; ii < m; ii++) {
                    relatedDepCandidates[ii].insert(i * m + ii);
                }
            }
            return relatedDepCandidates;
        }

        void testResults(ADFun<CGD>& fun, size_t m, size_t repeat, size_t n_loops) {
            /**
             * Generate operation graph
             */
            CodeHandler<double> h;
            size_t n2 = fun.Domain();

            std::vector<CGD> xx(n2);
            h.makeVariables(xx);
            for (size_t j = 0; j < n2; j++) {
                xx[j].setValue(j);
            }

            std::vector<CGD> yy = fun.Forward(0, xx);

            std::vector<std::set<size_t> > relatedDepCandidates = createRelatedDepCandidates(m, repeat);

            DependentPatternMatcher<double> matcher(relatedDepCandidates);

            std::vector<Loop<Base>*> loops = matcher.findLoops(yy, xx);
            std::cout << "loops: " << loops.size() << std::endl;
            ASSERT_EQ(loops.size(), n_loops);

            std::vector<EquationPattern<double>*> equations = matcher.getEquationPatterns();
            std::cout << "equation patterns: " << equations.size() << std::endl;
            ASSERT_EQ(equations.size(), m);

            std::auto_ptr<ADFun<CG<Base> > > newTape(matcher.createNewTape(yy, xx));

            // clean-up
            for (size_t l = 0; l < loops.size(); l++) {
                delete loops[l];
            }
        }

        void testSourceCodeGen(ADFun<CGD>& fun,
                               size_t m, size_t repeat,
                               const std::string& name,
                               JacobianADMode jacMode,
                               TEST_TYPE jacobian = MUST_PASS,
                               TEST_TYPE hessian = MUST_PASS,
                               bool reverseTwo = false) {
            std::vector<atomic_base<Base>*> atoms;
            testSourceCodeGen(fun, m, repeat, name, atoms, jacMode, jacobian, hessian, reverseTwo);
        }

        void testSourceCodeGen(ADFun<CGD>& fun,
                               size_t m, size_t repeat,
                               const std::string& name,
                               const std::vector<atomic_base<Base>*>& atoms,
                               JacobianADMode jacMode,
                               TEST_TYPE jacobian = MUST_PASS,
                               TEST_TYPE hessian = MUST_PASS,
                               bool reverseTwo = false) {

            std::string libBaseName = name;
            if (jacobian == MUST_PASS) {
                if (jacMode == FORWARD)libBaseName += "F";
                else if (jacMode == REVERSE)libBaseName += "R";
            }
            if (hessian == MUST_PASS && reverseTwo)
                libBaseName += "rev2";

            std::vector<std::set<size_t> > relatedDepCandidates = createRelatedDepCandidates(m, repeat);
            std::vector<double> xTypical(fun.Domain(), 0.9);
            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            CLangCompileModelHelper<double> compHelpL(fun, libBaseName + "DynamicWithLoops");
            compHelpL.setCreateForwardZero(true);
            compHelpL.setJacobianADMode(jacMode);
            compHelpL.setCreateJacobian(false);
            compHelpL.setCreateHessian(false);
            compHelpL.setCreateSparseJacobian(jacobian == MUST_PASS);
            compHelpL.setCreateSparseHessian(hessian == MUST_PASS);
            compHelpL.setCreateForwardOne(false);
            compHelpL.setCreateReverseOne(false);
            compHelpL.setCreateReverseTwo(reverseTwo);
            //compHelpL.setMaxAssignmentsPerFunc(maxAssignPerFunc);
            compHelpL.setRelatedDependents(relatedDepCandidates);
            compHelpL.setTypicalIndependentValues(xTypical);

            GccCompiler<double> compiler;
            compiler.setSourcesFolder("sources_" + libBaseName + "_1");

            CLangCompileDynamicHelper<double> compDynHelpL(compHelpL);
            std::auto_ptr<DynamicLib<double> > dynamicLibL(compDynHelpL.createDynamicLibrary(compiler));
            std::auto_ptr<DynamicLibModel<double> > modelL(dynamicLibL->model(libBaseName + "DynamicWithLoops"));
            for (size_t i = 0; i < atoms.size(); i++)
                modelL->addAtomicFunction(*atoms[i]);
            /**
             * Without the loops
             */
            CLangCompileModelHelper<double> compHelp(fun, libBaseName + "DynamicNoLoops");
            compHelp.setCreateForwardZero(true);
            compHelp.setJacobianADMode(jacMode);
            compHelp.setCreateJacobian(false);
            compHelp.setCreateHessian(false);
            compHelp.setCreateSparseJacobian(jacobian == MUST_PASS);
            compHelp.setCreateSparseHessian(hessian == MUST_PASS);
            compHelp.setCreateForwardOne(false);
            compHelp.setCreateReverseOne(false);
            compHelp.setCreateReverseTwo(reverseTwo);
            //compHelp.setMaxAssignmentsPerFunc(maxAssignPerFunc);

            compiler.setSourcesFolder("sources_" + libBaseName + "_1");

            CLangCompileDynamicHelper<double> compDynHelp(compHelp);
            compDynHelp.setLibraryName("modellibNoLoops");
            std::auto_ptr<DynamicLib<double> > dynamicLib(compDynHelp.createDynamicLibrary(compiler));

            /**
             * reference library
             */
            std::auto_ptr<DynamicLibModel<double> > model(dynamicLib->model(libBaseName + "DynamicNoLoops"));
            for (size_t i = 0; i < atoms.size(); i++)
                model->addAtomicFunction(*atoms[i]);

            /**
             * Compare results
             */
            size_t nFull = modelL->Domain();
            ASSERT_EQ(modelL->Domain(), model->Domain());
            ASSERT_EQ(modelL->Range(), model->Range());

            std::vector<double> x(nFull);
            for (size_t j = 0; j < nFull; j++) {
                x[j] = j + 1;
            }

            // test model (zero-order)
            if (compHelp.isCreateForwardZero()) {
                std::vector<double> yl = modelL->ForwardZero(x);
                std::vector<double> y = model->ForwardZero(x);
                compareValues(yl, y);
            }

            // test jacobian
            if (compHelp.isCreateSparseJacobian()) {
                compareVectorSetValues(modelL->JacobianSparsitySet(),
                                       model->JacobianSparsitySet());

                std::vector<double> jacl, jac;
                std::vector<size_t> rowsl, colsl, rows, cols;
                modelL->SparseJacobian(x, jacl, rowsl, colsl);
                model->SparseJacobian(x, jac, rows, cols);

                compareValues(jacl, jac);
            }

            // test hessian
            if (compHelp.isCreateSparseHessian()) {
                compareVectorSetValues(modelL->HessianSparsitySet(),
                                       model->HessianSparsitySet());

                std::vector<double> w(m * repeat);
                for (size_t i = 0; i < w.size(); i++) {
                    w[i] = 0.5 * (i + 1);
                }
                std::vector<double> hessl, hess;
                std::vector<size_t> rowsl, colsl, rows, cols;
                modelL->SparseHessian(x, w, hessl, rowsl, colsl);
                model->SparseHessian(x, w, hess, rows, cols);

                compareValues(hessl, hess);
            }


            if (jacobian == MUST_FAIL) {
                /** Make sure it fails */
                CLangCompileModelHelper<double> compHelpL(fun, libBaseName + "DynamicWithLoopsJacFail");
                compHelpL.setCreateForwardZero(false);
                compHelpL.setJacobianADMode(jacMode);
                compHelpL.setCreateSparseJacobian(true);
                compHelpL.setRelatedDependents(relatedDepCandidates);
                compHelpL.setTypicalIndependentValues(xTypical);

                GccCompiler<double> compiler;
                compiler.setSourcesFolder("sources_" + libBaseName + "_1");

                CLangCompileDynamicHelper<double> compDynHelpL(compHelpL);
                std::auto_ptr<DynamicLib<double> > dynamicLibL;
                ASSERT_THROW(dynamicLibL = std::auto_ptr<DynamicLib<double> >(compDynHelpL.createDynamicLibrary(compiler)), CGException);
            }

            if (hessian == MUST_FAIL) {
                /** Make sure it fails */
                CLangCompileModelHelper<double> compHelpL(fun, libBaseName + "DynamicWithLoopsHessFail");
                compHelpL.setCreateForwardZero(false);
                compHelpL.setCreateSparseHessian(true);
                compHelpL.setRelatedDependents(relatedDepCandidates);
                compHelpL.setTypicalIndependentValues(xTypical);

                GccCompiler<double> compiler;
                compiler.setSourcesFolder("sources_" + libBaseName + "_1");

                CLangCompileDynamicHelper<double> compDynHelpL(compHelpL);
                std::auto_ptr<DynamicLib<double> > dynamicLibL;
                ASSERT_THROW(dynamicLibL = std::auto_ptr<DynamicLib<double> >(compDynHelpL.createDynamicLibrary(compiler)), CGException);
            }
        }
    };
}

#endif