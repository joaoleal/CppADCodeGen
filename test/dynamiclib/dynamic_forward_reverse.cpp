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

#define _MODEL1

    class CppADCGDynamicForRevTest2 : public CppADCGTest {
    protected:
        const std::string _modelName;
        const static size_t n;
        const static size_t m;
        std::vector<double> x;
        ADFun<CGD>* _fun;
        DynamicLib<double>* _dynamicLib;
        DynamicLibModel<double>* _model;
    public:

        inline CppADCGDynamicForRevTest2(bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues),
            _modelName("model"),
            x(n),
            _fun(NULL),
            _dynamicLib(NULL),
            _model(NULL) {
        }

        virtual void SetUp() {
            // use a special object for source code generation
            typedef double Base;
            typedef CG<Base> CGD;
            typedef AD<CGD> ADCG;

            x[0] = 0.5;
            x[1] = 1.5;

            // independent variables
            std::vector<ADCG> u(n);
            for (size_t j = 0; j < n; j++)
                u[j] = x[j];

            CppAD::Independent(u);

            // dependent variable vector 
            std::vector<ADCG> Z(m);

            /**
             * create the CppAD tape as usual
             */
            Z[0] = 1.5 * x[0] + 1;
            Z[1] = 1.0 * x[1] + 2;

            // create f: U -> Z and vectors used for derivative calculations
            _fun = new ADFun<CGD>(u, Z);

            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            CLangCompileModelHelper<double> compHelp(*_fun, _modelName);

            compHelp.setCreateForwardZero(true);
            compHelp.setCreateForwardOne(true);
            compHelp.setCreateReverseOne(true);
            compHelp.setCreateReverseTwo(true);
            compHelp.setCreateSparseJacobian(true);
            compHelp.setCreateSparseHessian(true);

            GccCompiler<double> compiler;

            CLangCompileDynamicHelper<double> compDynHelp(compHelp);
            _dynamicLib = compDynHelp.createDynamicLibrary(compiler);
            _model = _dynamicLib->model(_modelName);

            // dimensions
            ASSERT_EQ(_model->Domain(), _fun->Domain());
            ASSERT_EQ(_model->Range(), _fun->Range());
        }

        virtual void TearDown() {
            delete _dynamicLib;
            _dynamicLib = NULL;
            delete _model;
            _model = NULL;
            delete _fun;
            _fun = NULL;
        }

    };

    /**
     * static data
     */
    const size_t CppADCGDynamicForRevTest2::n = 2;
    const size_t CppADCGDynamicForRevTest2::m = 2;

}

using namespace CppAD;
using namespace std;

TEST_F(CppADCGDynamicForRevTest2, SparseJacobian) {
    using namespace std;
    using std::vector;

    vector<CGD> xOrig(n);
    for (size_t i = 0; i < n; i++)
        xOrig[i] = x[i];

    const std::vector<bool> p = jacobianSparsity < std::vector<bool>, CGD > (*_fun);

    vector<CGD> jacOrig = _fun->SparseJacobian(xOrig, p);
    vector<double> jacCG = CppADCGDynamicForRevTest2::_model->SparseJacobian(x);

    compareValues(jacCG, jacOrig);
}

TEST_F(CppADCGDynamicForRevTest2, SparseHessian) {
    using namespace std;
    using std::vector;

    typedef double Base;
    typedef CppAD::CG<Base> CGD;

    vector<double> w(m);
    for (size_t i = 0; i < m; i++)
        w[i] = 1;

    vector<CGD> wOrig(m);
    for (size_t i = 0; i < m; i++)
        wOrig[i] = w[i];

    vector<CGD> xOrig(n);
    for (size_t i = 0; i < n; i++)
        xOrig[i] = x[i];

    vector<CGD> hessOrig = _fun->SparseHessian(xOrig, wOrig);
    vector<double> hessCG = CppADCGDynamicForRevTest2::_model->SparseHessian(x, w);

    compareValues(hessCG, hessOrig);
}