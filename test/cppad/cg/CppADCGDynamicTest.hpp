#ifndef CPPAD_CG_TEST_CPPADCGDYNAMICTEST_INCLUDED
#define CPPAD_CG_TEST_CPPADCGDYNAMICTEST_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
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
#include "CppADCGModelTest.hpp"
#include "gccCompilerFlags.hpp"

#ifdef CPPAD_CG_SYSTEM_LINUX
#include <dlfcn.h>
#endif

namespace CppAD {
namespace cg {

class CppADCGDynamicTest : public CppADCGModelTest {
public:
    using CGD = CG<double>;
    using ADCG = AD<CGD>;
protected:
    const std::string _name;
    bool _denseJacobian;
    bool _denseHessian;
    bool _forwardOne;
    bool _reverseOne;
    bool _reverseTwo;
    MultiThreadingType _multithread;
    bool _multithreadDisabled;
    ThreadPoolScheduleStrategy _multithreadScheduler;
public:

    explicit CppADCGDynamicTest(std::string testName,
                                bool verbose = false,
                                bool printValues = false) :
        CppADCGModelTest(verbose, printValues),
        _name(std::move(testName)),
        _denseJacobian(true),
        _denseHessian(true),
        _forwardOne(true),
        _reverseOne(true),
        _reverseTwo(true),
        _multithread(MultiThreadingType::NONE),
        _multithreadDisabled(false),
        _multithreadScheduler(ThreadPoolScheduleStrategy::DYNAMIC) {
    }

    virtual std::vector<ADCGD> model(const std::vector<ADCGD>& ind) = 0;

    virtual std::vector<ADCGD> model(const std::vector<ADCGD>& ind,
                                     const std::vector<ADCGD>& p) {
        return model(ind);
    }

    void testDynamicFull(std::vector<ADCG>& u,
                         const std::vector<double>& x,
                         size_t maxAssignPerFunc = 100,
                         double epsilonR = 1e-14,
                         double epsilonA = 1e-14) {
        std::vector<ADCG> p;
        const std::vector<double> xDynPar;

        testDynamicFull(u, p, x, xDynPar, maxAssignPerFunc, epsilonR, epsilonA);
    }

    void testDynamicFull(std::vector<ADCG>& u,
                         std::vector<ADCG>& p,
                         const std::vector<double>& x,
                         const std::vector<double>& xDynPar,
                         size_t maxAssignPerFunc = 100,
                         double epsilonR = 1e-14,
                         double epsilonA = 1e-14) {
        const std::vector<double> xNorm(x.size(), 1.0);
        const std::vector<double> eqNorm;

        testDynamicFull(u, p, x, xNorm, xDynPar, eqNorm, maxAssignPerFunc, epsilonR, epsilonA);
    }

    void testDynamicFull(std::vector<ADCG>& u,
                         const std::vector<double>& x,
                         const std::vector<double>& xNorm,
                         const std::vector<double>& eqNorm,
                         size_t maxAssignPerFunc = 100,
                         double epsilonR = 1e-14,
                         double epsilonA = 1e-14) {
        std::vector<ADCG> p;
        const std::vector<double> xDynPar;
        testDynamicFull(u, p, x, xNorm, xDynPar, eqNorm, maxAssignPerFunc, epsilonR, epsilonA);
    }

    void testDynamicFull(std::vector<ADCG>& ax,
                         std::vector<ADCG>& ap,
                         const std::vector<double>& x,
                         const std::vector<double>& xNorm,
                         const std::vector<double>& p,
                         const std::vector<double>& eqNorm,
                         size_t maxAssignPerFunc = 100,
                         double epsilonR = 1e-14,
                         double epsilonA = 1e-14) {
        ASSERT_EQ(ax.size(), x.size());
        ASSERT_EQ(x.size(), xNorm.size());

        using namespace std;

        // use a special object for source code generation
        // declare independent variables, dynamic parameters, starting recording
        size_t abort_op_index = 0;
        bool record_compare = true;
        CppAD::Independent(ax, abort_op_index, record_compare, ap);

        for (size_t i = 0; i < ax.size(); i++)
            ax[i] *= xNorm[i];

        // dependent variable vector
        std::vector<ADCG> ay = model(ax, ap);

        if (!eqNorm.empty()) {
            ASSERT_EQ(ay.size(), eqNorm.size());
            for (size_t i = 0; i < ay.size(); i++)
                ay[i] /= eqNorm[i];
        }

        /**
         * create the CppAD tape as usual
         */
        // create f: ax -> ay and vectors used for derivative calculations
        ADFun<CGD> fun;
        fun.Dependent(ay);

        /**
         * Create the dynamic library
         * (generate and compile source code)
         */
        ModelCSourceGen<double> compHelp(fun, _name + "dynamic");

        compHelp.setCreateForwardZero(true);
        compHelp.setCreateJacobian(_denseJacobian);
        compHelp.setCreateHessian(_denseHessian);
        compHelp.setCreateSparseJacobian(true);
        compHelp.setCreateSparseHessian(true);
        compHelp.setCreateForwardOne(_forwardOne);
        compHelp.setCreateReverseOne(_reverseOne);
        compHelp.setCreateReverseTwo(_reverseTwo);
        compHelp.setMaxAssignmentsPerFunc(maxAssignPerFunc);
        compHelp.setMultiThreading(true);

        ModelLibraryCSourceGen<double> compDynHelp(compHelp);
        compDynHelp.setMultiThreading(_multithread);

        SaveFilesModelLibraryProcessor<double>::saveLibrarySourcesTo(compDynHelp, "sources_" + _name + "_1");

        DynamicModelLibraryProcessor<double> libProc(compDynHelp);
        GccCompiler<double> compiler;
        //compiler.setSaveToDiskFirst(true); // useful to detect problem
        prepareTestCompilerFlags(compiler);
        if (compDynHelp.getMultiThreading() == MultiThreadingType::OPENMP) {
            compiler.addCompileFlag("-fopenmp");
            compiler.addCompileFlag("-pthread");
            compiler.addCompileLibFlag("-fopenmp");

#ifdef CPPAD_CG_SYSTEM_LINUX
            // this is required because the OpenMP implementation in GCC causes a segmentation fault on dlclose
            libProc.getOptions() ["dlOpenMode"] = std::to_string(RTLD_NOW | RTLD_NODELETE);
#endif
        } else if (compDynHelp.getMultiThreading() == MultiThreadingType::PTHREADS) {
            compiler.addCompileFlag("-pthread");
        }

        std::unique_ptr<DynamicLib<double>> dynamicLib = libProc.createDynamicLibrary(compiler);
        dynamicLib->setThreadPoolVerbose(this->verbose_);
        dynamicLib->setThreadNumber(2);
        dynamicLib->setThreadPoolDisabled(_multithreadDisabled);
        dynamicLib->setThreadPoolSchedulerStrategy(_multithreadScheduler);
        dynamicLib->setThreadPoolGuidedMaxWork(0.75);

        /**
         * test the library
         */
        std::unique_ptr<GenericModel<double>> model = dynamicLib->model(_name + "dynamic");
        ASSERT_TRUE(model != nullptr);

        testModelResults(*dynamicLib, *model, fun, x, p, epsilonR, epsilonA, _denseJacobian, _denseHessian);
    }

    void testDynamicCustomElements(std::vector<ADCG>& ax,
                                   const std::vector<double>& x,
                                   const std::vector<size_t>& jacRow, const std::vector<size_t>& jacCol,
                                   const std::vector<size_t>& hessRow, const std::vector<size_t>& hessCol) {
        std::vector<ADCG> ap;
        std::vector<double> p;

        testDynamicCustomElements(ax, ap, x, p, jacRow, jacCol, hessRow, hessCol);
    }

    void testDynamicCustomElements(std::vector<ADCG>& ax,
                                   std::vector<ADCG>& ap,
                                   const std::vector<double>& x,
                                   const std::vector<double>& p,
                                   const std::vector<size_t>& jacRow, const std::vector<size_t>& jacCol,
                                   const std::vector<size_t>& hessRow, const std::vector<size_t>& hessCol) {
        using namespace std;

        size_t abort_op_index = 0;
        bool record_compare = true;
        CppAD::Independent(ax, abort_op_index, record_compare, ap);

        // dependent variable vector
        std::vector<ADCG> ay = model(ax);

        /**
         * create the CppAD tape as usual
         */
        // create f: ax -> ay and vectors used for derivative calculations
        ADFun<CGD> fun(ax, ay);

        /**
         * Create the dynamic library
         * (generate and compile source code)
         */
        ModelCSourceGen<double> compHelp(fun, _name + "dynamic2");

        compHelp.setCreateForwardOne(_forwardOne);
        compHelp.setCreateReverseOne(_reverseOne);
        compHelp.setCreateReverseTwo(_reverseTwo);

        compHelp.setCreateSparseJacobian(true);
        compHelp.setCustomSparseJacobianElements(jacRow, jacCol);

        compHelp.setCreateSparseHessian(true);
        compHelp.setCustomSparseHessianElements(hessRow, hessCol);

        compHelp.setMultiThreading(true);

        ModelLibraryCSourceGen<double> compDynHelp(compHelp);
        compDynHelp.setMultiThreading(_multithread);

        SaveFilesModelLibraryProcessor<double>::saveLibrarySourcesTo(compDynHelp, "sources_" + _name + "_2");

        DynamicModelLibraryProcessor<double> dmlp(compDynHelp, "cppad_cg_model_2");

        GccCompiler<double> compiler;
        prepareTestCompilerFlags(compiler);
        if (compDynHelp.getMultiThreading() == MultiThreadingType::OPENMP) {
            compiler.addCompileFlag("-fopenmp");
            compiler.addCompileFlag("-pthread");
            compiler.addCompileLibFlag("-fopenmp");

#ifdef CPPAD_CG_SYSTEM_LINUX
            // this is required because the OpenMP implementation in GCC causes a segmentation fault on dlclose
            dmlp.getOptions() ["dlOpenMode"] = std::to_string(RTLD_NOW | RTLD_NODELETE);
#endif
        } else if (compDynHelp.getMultiThreading() == MultiThreadingType::PTHREADS) {
            compiler.addCompileFlag("-pthread");
        }

        std::unique_ptr<DynamicLib<double>> dynamicLib = dmlp.createDynamicLibrary(compiler);
        dynamicLib->setThreadPoolVerbose(this->verbose_);
        dynamicLib->setThreadNumber(2);
        dynamicLib->setThreadPoolDisabled(_multithreadDisabled);
        dynamicLib->setThreadPoolSchedulerStrategy(_multithreadScheduler);
        dynamicLib->setThreadPoolGuidedMaxWork(0.75);

        /**
         * test the library
         */
        std::unique_ptr<GenericModel<double>> model = dynamicLib->model(_name + "dynamic2");
        ASSERT_TRUE(model != nullptr);

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
        model->SparseJacobian(x, p, jacCGen, row, col);
        std::vector<CG<double> > jacSparse(row.size());

        std::vector<CGD> jac = fun.Jacobian(x2);
        for (size_t i = 0; i < row.size(); i++) {
            jacSparse[i] = jac[row[i] * x.size() + col[i]];
        }

        ASSERT_TRUE(compareValues(jacCGen, jacSparse));

        // sparse Hessian
        std::vector<double> w(ay.size(), 1.0);
        std::vector<double> hessCGen;
        model->SparseHessian(x, p, w, hessCGen, row, col);
        std::vector<CG<double> > hessSparse(row.size());

        std::vector<CGD> w2(ay.size(), 1.0);
        std::vector<CGD> hess = fun.Hessian(x2, w2);
        for (size_t i = 0; i < row.size(); i++) {
            hessSparse[i] = hess[row[i] * x.size() + col[i]];
        }

        ASSERT_TRUE(compareValues(hessCGen, hessSparse));
    }

};

} // END cg namespace
} // END CppAD namespace

#endif
