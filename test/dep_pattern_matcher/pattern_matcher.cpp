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

typedef double Base;
typedef CppAD::CG<Base> CGD;
typedef CppAD::AD<CGD> ADCGD;

namespace CppAD {

    class CppADCGPatternTest : public CppADCGTest {
    public:

        inline CppADCGPatternTest(bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues) {
        }

        void testPattern(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat),
                         size_t m,
                         size_t n,
                         size_t repeat,
                         size_t n_loops = 1,
                         bool createDynLib = false,
                         std::string name = "") {
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

            if (createDynLib) {
                testSourceCodeGen(fun, m, repeat, name, FORWARD);
                testSourceCodeGen(fun, m, repeat, name, REVERSE);
            } else {
                testResults(fun, m, repeat, n_loops);
            }
        }

        void testPatternWithAtomics(std::vector<ADCGD> (*model)(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<Base>*>& atoms),
                                    const std::vector<atomic_base<Base>* >& atoms,
                                    size_t m,
                                    size_t n,
                                    size_t repeat,
                                    size_t n_loops = 1,
                                    bool createDynLib = false,
                                    std::string name = "") {
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

            if (createDynLib) {
                testSourceCodeGen(fun, m, repeat, name, atoms, FORWARD);
                testSourceCodeGen(fun, m, repeat, name, atoms, REVERSE);
            } else {
                testResults(fun, m, repeat, n_loops);
            }

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
                               JacobianADMode jacMode) {
            std::vector<atomic_base<Base>*> atoms;
            testSourceCodeGen(fun, m, repeat, name, atoms, jacMode);
        }

        void testSourceCodeGen(ADFun<CGD>& fun,
                               size_t m, size_t repeat,
                               const std::string& name,
                               const std::vector<atomic_base<Base>*>& atoms,
                               JacobianADMode jacMode) {

            std::string libBaseName = name;
            if (jacMode == FORWARD)libBaseName += "F";
            else if (jacMode == REVERSE)libBaseName += "R";

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
            compHelpL.setCreateSparseJacobian(true);
            compHelpL.setCreateSparseHessian(false);
            compHelpL.setCreateForwardOne(false);
            compHelpL.setCreateReverseOne(false);
            compHelpL.setCreateReverseTwo(false);
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
            compHelp.setCreateSparseJacobian(true);
            compHelp.setCreateSparseHessian(false);
            compHelp.setCreateForwardOne(false);
            compHelp.setCreateReverseOne(false);
            compHelp.setCreateReverseTwo(false);
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

            if (compHelp.isCreateForwardZero()) {
                std::vector<double> yl = modelL->ForwardZero(x);
                std::vector<double> y = model->ForwardZero(x);
                compareValues(yl, y);
            }

            if (compHelp.isCreateSparseJacobian()) {
                compareVectorSetValues(modelL->JacobianSparsitySet(),
                                       model->JacobianSparsitySet());

                std::vector<double> jacl, jac;
                std::vector<size_t> rowsl, colsl, rows, cols;
                modelL->SparseJacobian(x, jacl, rowsl, colsl);
                model->SparseJacobian(x, jac, rows, cols);

                compareValues(jacl, jac);
            }

        }
    };
}

using namespace CppAD;

std::vector<ADCGD> modelCommonTmp2(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp1 = x[1] * sin(x[1] + 4);
    ADCGD tmp2 = 2.5 * x[1];
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = tmp1 * cos(x[i * m]) / tmp2;
        y[i * m + 1] = x[i * m + 1] * x[i * m];
    }

    return y;
}

TEST_F(CppADCGPatternTest, CommonTmp2) {
    using namespace CppAD;
    size_t m = 2;
    size_t n = 2;

    testPattern(modelCommonTmp2, m, n, 6);

    testPattern(modelCommonTmp2, m, n, 6, 1, true, "modelCommonTmp2");
}

std::vector<ADCGD> model0(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * m]);
        y[i * m + 1] = x[i * m + 1] * x[i * m];
    }

    return y;
}

TEST_F(CppADCGPatternTest, DependentPatternMatcherDetached) {
    using namespace CppAD;
    size_t m = 2;
    size_t n = 2;

    testPattern(model0, m, n, 6);

    testPattern(model0, m, n, 6, 1, true, "model0");
}

std::vector<ADCGD> model1(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * m]) + x[1] * x[2];
        y[i * m + 1] = x[i * m + 1] * x[i * m];
    }

    return y;
}

TEST_F(CppADCGPatternTest, DependentPatternMatcher) {
    using namespace CppAD;
    size_t m = 2;
    size_t n = 2;

    testPattern(model1, m, n, 6);

    testPattern(model1, m, n, 6, 1, true, "model1");
}

std::vector<ADCGD> model4Eq(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 4;
    size_t m2 = repeat * m;

    assert(x.size() == m2);

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * m]) + x[1] * log(x[2]);
        y[i * m + 1] = x[i * m + 1] * x[i * m];
        y[i * m + 2] = x[i * m + 1] * x[i * m + 2];
        y[i * m + 3] = 5;
    }

    return y;
}

TEST_F(CppADCGPatternTest, Matcher4Eq) {
    using namespace CppAD;
    size_t m = 4;
    size_t n = 4;

    testPattern(model4Eq, m, n, 6);

    testPattern(model4Eq, m, n, 6, 1, true, "model4Eq");
}

std::vector<ADCGD> modelCommonTmp(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp = x[1] * x[2];
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * m]) + tmp;
        y[i * m + 1] = x[i * m + 1] * x[i * m] + tmp;
    }

    return y;
}

TEST_F(CppADCGPatternTest, CommonTmp) {
    using namespace CppAD;
    size_t m = 2;
    size_t n = 2;

    testPattern(modelCommonTmp, m, n, 6);

    testPattern(modelCommonTmp, m, n, 6, 1, true, "modelCommonTmp");
}

std::vector<ADCGD> model4(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        ADCGD tmp = x[1] * x[i];
        y[i * m] = cos(x[i * m]) + tmp;
        y[i * m + 1] = x[i * m + 1] * x[i * m] + tmp;
    }

    return y;
}

TEST_F(CppADCGPatternTest, IndexedTmp) {
    using namespace CppAD;

    size_t m = 2;
    size_t n = 2;

    testPattern(model4, m, n, 6);

    testPattern(model4, m, n, 6, 1, true, "indexedTmp");
}

std::vector<ADCGD> model5(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        ADCGD tmp = x[1] * x[i];
        y[i * m] = cos(x[i * m]) + tmp;

        if (i == 1) {
            for (size_t i2 = 0; i2 < repeat; i2++) {
                y[i2 * m + 1] = x[i2 * m + 1] * x[i2 * m] + tmp;
            }
        }
    }

    return y;
}

TEST_F(CppADCGPatternTest, DependentPatternMatcher5) {
    using namespace CppAD;

    size_t m = 2;
    size_t n = 2;

    testPattern(model5, m, n, 6, 2);

    testPattern(model5, m, n, 6, 2, true, "model5");
}

std::vector<ADCGD> modelAtomic(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<double>*>& atoms) {
    size_t m = 2;
    size_t m2 = repeat * m;

    CGAbstractAtomicFun<Base>& atomic0 = *atoms[0];

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    std::vector<ADCGD> ax(2), ay(1);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * m]);

        ax[0] = x[i * m];
        ax[1] = x[i * m + 1];
        atomic0(ax, ay);
        y[i * m + 1] = ay[0];
    }

    return y;
}

void atomicFunction(const std::vector<AD<double> >& x, std::vector<AD<double> >& y) {
    y[0] = x[1] * x[0];
}

TEST_F(CppADCGPatternTest, Atomic) {
    using namespace CppAD;

    size_t m = 2;
    size_t n = 2;

    // create atomic function
    std::vector<AD<double> > y(1);
    std::vector<AD<double> > x(2);

    checkpoint<double> atomicfun("atomicFunc", atomicFunction, x, y);
    std::vector<atomic_base<double>*> atomics(1);
    atomics[0] = &atomicfun;


    testPatternWithAtomics(modelAtomic, atomics, m, n, 6);

    testPatternWithAtomics(modelAtomic, atomics, m, n, 6, 1, true, "modelAtomic");
}