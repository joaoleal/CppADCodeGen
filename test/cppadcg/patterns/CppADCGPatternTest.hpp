#ifndef CPPAD_CG_TEST_CPPADCGPATTERNTEST_INCLUDED
#define	CPPAD_CG_TEST_CPPADCGPATTERNTEST_INCLUDED
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

    class CppADCGPatternTest : public CppADCGTest {
    public:
        typedef double Base;
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCGD;
    protected:
        bool testZeroOrder_;
        bool testJacobian_;
        bool testHessian_;
        std::vector<Base> xNorm_;
        std::vector<Base> eqNorm_;
        Base epsilonA_;
        Base epsilonR_;
        Base hessianEpsilonA_;
        Base hessianEpsilonR_;
        std::vector<std::set<size_t> > customJacSparsity_;
        std::vector<std::set<size_t> > customHessSparsity_;
    public:

        inline CppADCGPatternTest(bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues),
            testZeroOrder_(true),
            testJacobian_(true),
            testHessian_(true),
            epsilonA_(std::numeric_limits<Base>::epsilon() * 1e2),
            epsilonR_(std::numeric_limits<Base>::epsilon() * 1e2),
            hessianEpsilonA_(std::numeric_limits<Base>::epsilon() * 1e2),
            hessianEpsilonR_(std::numeric_limits<Base>::epsilon() * 1e2) {
            //this->verbose_ = true;
        }

        void testPatternDetection(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat),
                                  size_t m,
                                  size_t n,
                                  size_t repeat,
                                  size_t n_loops = 1) {
            std::vector<std::vector<std::set<size_t> > > loops(n_loops);
            testPatternDetection(model, m, n, repeat, loops);
        }

        void testPatternDetection(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat),
                                  size_t m,
                                  size_t n,
                                  size_t repeat,
                                  const std::vector<std::vector<std::set<size_t> > >& loops) {
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

            testPatternDetectionResults(fun, m, repeat, loops);
        }

        void testLibCreation(const std::string& libName,
                             std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat),
                             size_t m,
                             size_t n,
                             size_t repeat,
                             size_t mExtra = 0) {
            size_t n2 = repeat * n;
            std::vector<Base> x(n2);
            for (size_t j = 0; j < n2; j++)
                x[j] = 0.5 * (j + 1);

            testLibCreation(libName, model, m, repeat, mExtra, x);
        }

        void testLibCreation(const std::string& libName,
                             std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat),
                             size_t m,
                             size_t repeat,
                             size_t mExtra,
                             const std::vector<Base>& xb) {
            using namespace CppAD;

            /**
             * Tape model
             */
            std::vector<ADCGD> x(xb.size());
            for (size_t j = 0; j < xb.size(); j++)
                x[j] = xb[j];
            CppAD::Independent(x);
            if (xNorm_.size() > 0) {
                ASSERT_EQ(x.size(), xNorm_.size());
                for (size_t j = 0; j < x.size(); j++)
                    x[j] *= xNorm_[j];
            }

            std::vector<ADCGD> y = (*model)(x, repeat);
            if (eqNorm_.size() > 0) {
                ASSERT_EQ(y.size(), eqNorm_.size());
                for (size_t i = 0; i < y.size(); i++)
                    y[i] /= eqNorm_[i];
            }

            ADFun<CGD> fun;
            fun.Dependent(y);

            testSourceCodeGen(fun, m, repeat, mExtra, libName, xb, FORWARD, testJacobian_, testHessian_);
            if (testJacobian_) {
                testSourceCodeGen(fun, m, repeat, mExtra, libName, xb, FORWARD, true, false, true);
                testSourceCodeGen(fun, m, repeat, mExtra, libName, xb, REVERSE, true, false);
                testSourceCodeGen(fun, m, repeat, mExtra, libName, xb, REVERSE, true, false, true);
            }

            if (testHessian_) {
                testSourceCodeGen(fun, m, repeat, mExtra, libName, xb, FORWARD, false, true, false, true);
            }

        }

        void testPatternDetectionWithAtomics(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<Base>*>& atoms),
                                             const std::vector<atomic_base<Base>* >& atoms,
                                             size_t m,
                                             size_t n,
                                             size_t repeat,
                                             size_t n_loops = 1) {
            size_t n2 = repeat * n;
            std::vector<Base> x(n2);
            for (size_t j = 0; j < n2; j++)
                x[j] = 0.5;

            testPatternDetectionWithAtomics(model, atoms, m, x, repeat, n_loops);
        }

        void testPatternDetectionWithAtomics(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<Base>*>& atoms),
                                             const std::vector<atomic_base<Base>* >& atoms,
                                             size_t m,
                                             const std::vector<Base>& xb,
                                             size_t repeat,
                                             size_t n_loops = 1) {

            std::vector<std::vector<std::set<size_t> > > loops(n_loops);
            testPatternDetectionWithAtomics(model, atoms, m, xb, repeat, loops);
        }

        void testPatternDetectionWithAtomics(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<Base>*>& atoms),
                                             const std::vector<atomic_base<Base>* >& atoms,
                                             size_t m,
                                             const std::vector<Base>& xb,
                                             size_t repeat,
                                             const std::vector<std::vector<std::set<size_t> > >& loops) {

            using namespace CppAD;

            std::vector<CGAbstractAtomicFun<double>*> atomics(atoms.size());
            for (size_t a = 0; a < atoms.size(); a++) {
                atomics[a] = new CGAtomicFun<Base>(*atoms[a], true);
            }

            //size_t m2 = repeat * m;
            size_t n2 = xb.size();

            /**
             * Tape model
             */
            std::vector<ADCGD> x(n2);
            for (size_t j = 0; j < n2; j++)
                x[j] = xb[j];
            CppAD::Independent(x);
            if (xNorm_.size() > 0) {
                ASSERT_EQ(x.size(), xNorm_.size());
                for (size_t j = 0; j < x.size(); j++)
                    x[j] *= xNorm_[j];
            }

            std::vector<ADCGD> y = (*model)(x, repeat, atomics);
            if (eqNorm_.size() > 0) {
                ASSERT_EQ(y.size(), eqNorm_.size());
                for (size_t i = 0; i < y.size(); i++)
                    y[i] /= eqNorm_[i];
            }

            ADFun<CGD> fun;
            fun.Dependent(y);

            testPatternDetectionResults(fun, m, repeat, loops);

            for (size_t a = 0; a < atomics.size(); a++) {
                delete atomics[a];
            }
        }

        void testLibCreationWithAtomics(const std::string& name,
                                        std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<Base>*>& atoms),
                                        const std::vector<atomic_base<Base>* >& atoms,
                                        size_t m,
                                        size_t n,
                                        size_t repeat) {
            size_t n2 = repeat * n;
            std::vector<Base> x(n2);
            for (size_t j = 0; j < n2; j++)
                x[j] = 0.5;

            testLibCreationWithAtomics(name, model, atoms, m, x, repeat);
        }

        void testLibCreationWithAtomics(const std::string& name,
                                        std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<Base>*>& atoms),
                                        const std::vector<atomic_base<Base>* >& atoms,
                                        size_t m,
                                        const std::vector<Base>& xb,
                                        size_t repeat) {
            using namespace CppAD;

            SmartVectorPointer<CGAbstractAtomicFun<double> > atomics(atoms.size());
            for (size_t a = 0; a < atoms.size(); a++) {
                atomics.v[a] = new CGAtomicFun<Base>(*atoms[a], true);
            }

            //size_t m2 = repeat * m;
            size_t n2 = xb.size();

            /**
             * Tape model
             */
            std::vector<ADCGD> x(n2);
            for (size_t j = 0; j < n2; j++)
                x[j] = xb[j];
            CppAD::Independent(x);

            std::vector<ADCGD> y = (*model)(x, repeat, atomics.v);

            ADFun<CGD> fun;
            fun.Dependent(y);

            size_t mExtra = 0;
            testSourceCodeGen(fun, m, repeat, mExtra, name, atoms, xb, FORWARD, testJacobian_, testHessian_);
            if (testJacobian_) {
                testSourceCodeGen(fun, m, repeat, mExtra, name, atoms, xb, FORWARD, true, false, true);
                testSourceCodeGen(fun, m, repeat, mExtra, name, atoms, xb, REVERSE, true, false);
                testSourceCodeGen(fun, m, repeat, mExtra, name, atoms, xb, REVERSE, true, false, true);
            }

            if (testHessian_) {
                testSourceCodeGen(fun, m, repeat, mExtra, name, atoms, xb, FORWARD, false, true, false, true);
            }
        }

        std::vector<std::set<size_t> > createRelatedDepCandidates(size_t m,
                                                                  size_t repeat) {
            std::vector<std::set<size_t> > relatedDepCandidates(m);
            for (size_t i = 0; i < repeat; i++) {
                for (size_t ii = 0; ii < m; ii++) {
                    relatedDepCandidates[ii].insert(i * m + ii);
                }
            }
            return relatedDepCandidates;
        }

        void testPatternDetectionResults(ADFun<CGD>& fun,
                                         size_t m,
                                         size_t repeat,
                                         const std::vector<std::vector<std::set<size_t> > >& loops) {
            using namespace std;

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

            DependentPatternMatcher<double> matcher(relatedDepCandidates, yy, xx);

            LoopFreeModel<Base>* nonLoopTape;
            SmartSetPointer<LoopModel<Base> > loopTapes;
            matcher.generateTapes(nonLoopTape, loopTapes.s);

            delete nonLoopTape;

            //std::cout << "loops: " << matcher.getLoops().size() << std::endl;
            ASSERT_EQ(loopTapes.s.size(), loops.size());

            //std::cout << "equation patterns: " << matcher.getEquationPatterns().size() << std::endl;

            /**
             * order loops and equation patterns by the lowest used dependent
             */
            //  - calculated
            map<size_t, map<size_t, set<size_t> > > orderedCalcLoops;

            const std::vector<Loop<Base>*>& calcLoops = matcher.getLoops();
            for (size_t l = 0; l < calcLoops.size(); l++) {
                Loop<Base>* loop = calcLoops[l];
                size_t minDep = std::numeric_limits<size_t>::max();

                map<size_t, set<size_t> > dependents;
                set<EquationPattern<Base>*>::const_iterator iteq;
                for (iteq = loop->equations.begin(); iteq != loop->equations.end(); ++iteq) {
                    EquationPattern<Base>* eq = *iteq;
                    size_t minEqDep = *eq->dependents.begin();
                    minDep = std::min(minDep, *eq->dependents.begin());
                    dependents[minEqDep] = eq->dependents;
                }

                orderedCalcLoops[minDep] = dependents;
            }

            //  - expected
            bool defined = false;
            map<size_t, map<size_t, set<size_t> > > orderedExpectedLoops;
            for (size_t l = 0; l < loops.size(); l++) {
                const std::vector<set<size_t> >& eqPatterns = loops[l];

                if (!eqPatterns.empty()) {
                    defined = true;
                    /**
                     * check every equation
                     */
                    size_t minDep = std::numeric_limits<size_t>::max();

                    map<size_t, set<size_t> > dependents;
                    for (size_t eq = 0; eq < eqPatterns.size(); eq++) {
                        size_t minEqDep = *eqPatterns[eq].begin();
                        minDep = std::min(minDep, *eqPatterns[eq].begin());
                        dependents[minEqDep] = eqPatterns[eq];
                    }
                    orderedExpectedLoops[minDep] = dependents;
                }
            }

            if (defined) {
                map<size_t, map<size_t, set<size_t> > >::const_iterator itLexp = orderedExpectedLoops.begin();
                map<size_t, map<size_t, set<size_t> > >::const_iterator itLcalc = orderedCalcLoops.begin();
                for (; itLexp != orderedExpectedLoops.end(); ++itLexp, ++itLcalc) {
                    ASSERT_EQ(itLexp->first, itLcalc->first);
                    ASSERT_EQ(itLexp->second.size(), itLcalc->second.size());

                    map<size_t, set<size_t> >::const_iterator itEexp = itLexp->second.begin();
                    map<size_t, set<size_t> >::const_iterator itEcalc = itLcalc->second.begin();
                    for (; itEexp != itLexp->second.end(); ++itEexp, ++itEcalc) {
                        ASSERT_EQ(itEexp->first, itEcalc->first);
                        ASSERT_TRUE(itEexp->second == itEcalc->second);
                    }
                }

            } else {
                ASSERT_EQ(matcher.getEquationPatterns().size(), m);
                for (size_t eq = 0; eq < matcher.getEquationPatterns().size(); eq++) {
                    EquationPattern<Base>* eqp = matcher.getEquationPatterns()[eq];
                    ASSERT_EQ(eqp->dependents.size(), repeat);
                }
            }
        }

        void testSourceCodeGen(ADFun<CGD>& fun,
                               size_t m, size_t repeat,
                               size_t mExtra,
                               const std::string& name,
                               const std::vector<Base>& xTypical,
                               JacobianADMode jacMode,
                               bool jacobian = true,
                               bool hessian = true,
                               bool forReverseOne = false,
                               bool reverseTwo = false) {
            std::vector<atomic_base<Base>*> atoms;
            testSourceCodeGen(fun, m, repeat, mExtra, name, atoms, xTypical, jacMode, jacobian, hessian, forReverseOne, reverseTwo);
        }

        void testSourceCodeGen(ADFun<CGD>& fun,
                               size_t m, size_t repeat,
                               size_t mExtra,
                               const std::string& name,
                               const std::vector<atomic_base<Base>*>& atoms,
                               const std::vector<Base>& xTypical,
                               JacobianADMode jacMode,
                               bool jacobian = true,
                               bool hessian = true,
                               bool forReverseOne = false,
                               bool reverseTwo = false) {

            std::vector<std::set<size_t> > relatedDepCandidates = createRelatedDepCandidates(m, repeat);

            testSourceCodeGen(fun, relatedDepCandidates, name, atoms, xTypical,
                              jacMode, jacobian, hessian, forReverseOne, reverseTwo);
        }

        void testSourceCodeGen(ADFun<CGD>& fun,
                               const std::vector<std::set<size_t> >& relatedDepCandidates,
                               const std::string& name,
                               const std::vector<atomic_base<Base>*>& atoms,
                               const std::vector<Base>& xTypical,
                               JacobianADMode jacMode,
                               bool jacobian = true,
                               bool hessian = true,
                               bool forReverseOne = false,
                               bool reverseTwo = false) {

            bool loadModels = this->testZeroOrder_ || jacobian || hessian;

            std::string libBaseName = name;
            if (jacobian) {
                if (!forReverseOne) libBaseName += "d";
                if (jacMode == FORWARD) libBaseName += "F";
                else if (jacMode == REVERSE) libBaseName += "R";
            }
            if (hessian && reverseTwo)
                libBaseName += "rev2";


            assert(fun.Domain() == xTypical.size());
            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            CLangCompileModelHelper<double> compHelpL(fun, libBaseName + "Loops");
            compHelpL.setCreateForwardZero(true);
            compHelpL.setJacobianADMode(jacMode);
            compHelpL.setCreateJacobian(false);
            compHelpL.setCreateHessian(false);
            compHelpL.setCreateSparseJacobian(jacobian);
            compHelpL.setCreateSparseHessian(hessian);
            compHelpL.setCreateForwardOne(forReverseOne && jacMode == FORWARD);
            compHelpL.setCreateReverseOne(forReverseOne && jacMode == REVERSE);
            compHelpL.setCreateReverseTwo(reverseTwo);
            //compHelpL.setMaxAssignmentsPerFunc(maxAssignPerFunc);
            compHelpL.setRelatedDependents(relatedDepCandidates);
            compHelpL.setTypicalIndependentValues(xTypical);
            compHelpL.setParameterPrecision(std::numeric_limits<Base>::digits10 + 4);

            if (!customJacSparsity_.empty())
                compHelpL.setCustomSparseJacobianElements(customJacSparsity_);

            if (!customHessSparsity_.empty())
                compHelpL.setCustomSparseHessianElements(customHessSparsity_);

            GccCompiler<double> compiler;
            std::vector<std::string> flags;
            flags.push_back("-O0");
            flags.push_back("-g");
            flags.push_back("-ggdb");
            flags.push_back("-D_FORTIFY_SOURCE=2");
            compiler.setCompileFlags(flags);
            compiler.setSourcesFolder("sources_" + libBaseName);

            CLangCompileDynamicHelper<double> compDynHelpL(compHelpL);
            compDynHelpL.setVerbose(this->verbose_);
            std::auto_ptr<DynamicLib<double> > dynamicLibL(compDynHelpL.createDynamicLibrary(compiler));
            std::auto_ptr<DynamicLibModel<double> > modelL;
            if (loadModels) {
                modelL.reset(dynamicLibL->model(libBaseName + "Loops"));
                ASSERT_TRUE(modelL.get() != NULL);
                for (size_t i = 0; i < atoms.size(); i++)
                    modelL->addAtomicFunction(*atoms[i]);
            }
            /**
             * Without the loops
             */
            CLangCompileModelHelper<double> compHelp(fun, libBaseName + "NoLoops");
            compHelp.setCreateForwardZero(testZeroOrder_);
            compHelp.setJacobianADMode(jacMode);
            compHelp.setCreateJacobian(false);
            compHelp.setCreateHessian(false);
            compHelp.setCreateSparseJacobian(jacobian);
            compHelp.setCreateSparseHessian(hessian);
            compHelp.setCreateForwardOne(false);
            compHelp.setCreateReverseOne(false);
            compHelp.setCreateReverseTwo(reverseTwo);
            compHelp.setTypicalIndependentValues(xTypical);
            compHelp.setParameterPrecision(std::numeric_limits<Base>::digits10 + 4);
            //compHelp.setMaxAssignmentsPerFunc(maxAssignPerFunc);

            if (!customJacSparsity_.empty())
                compHelp.setCustomSparseJacobianElements(customJacSparsity_);

            if (!customHessSparsity_.empty())
                compHelp.setCustomSparseHessianElements(customHessSparsity_);

            compiler.setSourcesFolder("sources_" + libBaseName);

            CLangCompileDynamicHelper<double> compDynHelp(compHelp);
            compDynHelp.setVerbose(this->verbose_);
            compDynHelp.setLibraryName("modelLibNoLoops");
            std::auto_ptr<DynamicLib<double> > dynamicLib(compDynHelp.createDynamicLibrary(compiler));

            /**
             * reference library
             */
            std::auto_ptr<DynamicLibModel<double> > model;
            if (loadModels) {
                model.reset(dynamicLib->model(libBaseName + "NoLoops"));
                for (size_t i = 0; i < atoms.size(); i++)
                    model->addAtomicFunction(*atoms[i]);
            }

            if (!loadModels)
                return;

            /**
             * Compare results
             */
            ASSERT_EQ(modelL->Domain(), model->Domain());
            ASSERT_EQ(modelL->Range(), model->Range());

            std::vector<double> x = xTypical;

            // test model (zero-order)
            if (compHelp.isCreateForwardZero()) {
                std::vector<double> yl = modelL->ForwardZero(x);
                std::vector<double> y = model->ForwardZero(x);
                ASSERT_TRUE(compareValues(yl, y, epsilonR_, epsilonA_));
            }

            // test Jacobian
            if (compHelp.isCreateSparseJacobian()) {
                compareVectorSetValues(modelL->JacobianSparsitySet(),
                                       model->JacobianSparsitySet());

                std::vector<double> jacl, jac;
                std::vector<size_t> rowsl, colsl, rows, cols;
                modelL->SparseJacobian(x, jacl, rowsl, colsl);
                model->SparseJacobian(x, jac, rows, cols);

                ASSERT_TRUE(compareValues(jacl, jac, epsilonR_, epsilonA_));
            }

            // test Hessian
            if (compHelp.isCreateSparseHessian()) {
                compareVectorSetValues(modelL->HessianSparsitySet(),
                                       model->HessianSparsitySet());

                std::vector<double> w(fun.Range());
                for (size_t i = 0; i < w.size(); i++) {
                    w[i] = 0.5 * (i + 1);
                }
                std::vector<double> hessl, hess;
                std::vector<size_t> rowsl, colsl, rows, cols;
                modelL->SparseHessian(x, w, hessl, rowsl, colsl);
                model->SparseHessian(x, w, hess, rows, cols);

                ASSERT_TRUE(compareValues(hessl, hess, hessianEpsilonR_, hessianEpsilonA_));
            }

        }
    };
}

#endif