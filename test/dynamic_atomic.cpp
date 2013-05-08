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

#define MODEL1

    class CppADCGDynamicAtomicTest : public CppADCGTest {
    protected:
        const static std::string MODEL_NAME;
        const static size_t n;
        const static size_t m;
        CppAD::vector<double> x;
        ADFun<CGD>* _fun;
        DynamicLib<double>* _dynamicLib;
    public:
        static DynamicLibModel<double>* MODEL;
    public:

        inline CppADCGDynamicAtomicTest(bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues),
            x(n),
            _dynamicLib(NULL) {
        }

        virtual void SetUp() {
            // use a special object for source code generation
            typedef double Base;
            typedef CG<Base> CGD;
            typedef AD<CGD> ADCG;

#ifdef MODEL1
            for (size_t j = 0; j < n; j++)
                x[j] = j + 2;
#else
            x[0] = 0.5;
#endif
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
#ifdef MODEL1
            Z[0] = cos(u[0]);
            Z[1] = u[1] * u[2] + sin(u[0]);
            Z[2] = u[2] * u[2] + sin(u[1]);
            Z[3] = u[0] / u[2] + u[1] * u[2] + 5.0;
#else
            Z[0] = 1.0 / u[0];
#endif
            // create f: U -> Z and vectors used for derivative calculations
            _fun = new ADFun<CGD>(u, Z);

            /**
             * Create the dynamic library
             * (generate and compile source code)
             */
            CLangCompileModelHelper<double> compHelp(*_fun, MODEL_NAME);

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
            MODEL = _dynamicLib->model(MODEL_NAME);

            // dimensions
            ASSERT_EQ(MODEL->Domain(), _fun->Domain());
            ASSERT_EQ(MODEL->Range(), _fun->Range());
        }

        virtual void TearDown() {
            delete _dynamicLib;
            _dynamicLib = NULL;
            delete MODEL;
            MODEL = NULL;
            delete _fun;
            _fun = NULL;
        }

    };

    /**
     * static data
     */
    const std::string CppADCGDynamicAtomicTest::MODEL_NAME = "dynamicAtomic";
#ifdef MODEL1
    const size_t CppADCGDynamicAtomicTest::n = 3;
    const size_t CppADCGDynamicAtomicTest::m = 4;
#else
    const size_t CppADCGDynamicAtomicTest::n = 1;
    const size_t CppADCGDynamicAtomicTest::m = 1;
#endif
    DynamicLibModel<double>* CppADCGDynamicAtomicTest::MODEL = NULL;

    /**
     * User Defined Atomic AD Functions
     */
    bool atomic_use_lib_forward(size_t id,
                                size_t k,
                                size_t n,
                                size_t m,
                                const vector<bool>& vx,
                                vector<bool>& vy,
                                const vector<double>& tx,
                                vector<double>& ty) {

        if (k == 0) {
            CppADCGDynamicAtomicTest::MODEL->ForwardZero(vx, vy, tx, ty);
            return true;
        } else if (k == 1) {
            CppADCGDynamicAtomicTest::MODEL->ForwardOne(tx, ty);
            return true;
        }

        return false;
    }

    bool atomic_use_lib_reverse(size_t id,
                                size_t k,
                                size_t n,
                                size_t m,
                                const vector<double>& tx,
                                const vector<double>& ty,
                                vector<double>& px,
                                const vector<double>& py) {
        if (k == 0) {
            CppADCGDynamicAtomicTest::MODEL->ReverseOne(tx, ty, px, py);
            return true;
        } else if (k == 1) {
            CppADCGDynamicAtomicTest::MODEL->ReverseTwo(tx, ty, px, py);
            return true;
        }

        return false;
    }

    bool atomic_use_lib_for_jac_sparse(size_t id,
                                       size_t n,
                                       size_t m,
                                       size_t q,
                                       const CppAD::vector< std::set<size_t> >& r,
                                       CppAD::vector< std::set<size_t> >& s) {
        const std::vector<std::set<size_t> > jacSparsity = CppADCGDynamicAtomicTest::MODEL->JacobianSparsitySet();
        multMatrixMatrixSparsity(jacSparsity, r, s, q);
        return true;
    }

    bool atomic_use_lib_rev_jac_sparse(size_t id,
                                       size_t n,
                                       size_t m,
                                       size_t q,
                                       CppAD::vector< std::set<size_t> >& r,
                                       const CppAD::vector< std::set<size_t> >& s) {
        const std::vector<std::set<size_t> > jacSparsity = CppADCGDynamicAtomicTest::MODEL->JacobianSparsitySet();
        multMatrixMatrixSparsity(s, jacSparsity, r, n);
        return true;
    }

    bool atomic_use_lib_rev_hes_sparse(size_t id,
                                       size_t n,
                                       size_t m,
                                       size_t q,
                                       const CppAD::vector< std::set<size_t> >& r,
                                       const CppAD::vector<bool>& s,
                                       CppAD::vector<bool>& t,
                                       const CppAD::vector< std::set<size_t> >& u,
                                       CppAD::vector< std::set<size_t> >& v) {
        const std::vector<std::set<size_t> > jacSparsity = CppADCGDynamicAtomicTest::MODEL->JacobianSparsitySet();

        /**
         *  V(x)  =  f^1^T(x) U(x)  +  Sum(  s(x)i  f^2(x)  R(x)   )
         */
        // f^1^T(x) U(x)
        multMatrixTransMatrixSparsity(jacSparsity, u, v, q);

        // Sum(  s(x)i  f^2(x)  R(x)   )
        bool allSelected = true;
        for (size_t i = 0; i < s.size(); i++) {
            if (!s[i]) {
                allSelected = false;
                break;
            }
        }

        std::vector<std::set<size_t> > sparsitySF2R(n);
        if (allSelected) {
            sparsitySF2R = CppADCGDynamicAtomicTest::MODEL->HessianSparsitySet();
        } else {
            for (size_t i = 0; i < s.size(); i++) {
                if (s[i]) {
                    const std::vector<std::set<size_t> > sparsity = CppADCGDynamicAtomicTest::MODEL->HessianSparsitySet(i);
                    addMatrixSparsity(sparsity, sparsitySF2R);
                }
            }
        }

        multMatrixMatrixSparsity(sparsitySF2R, r, v, q);

        /**
         * S(x) * f^1(x)
         */
        for (size_t i = 0; i < s.size(); i++) {
            if (s[i]) {
                std::set<size_t>::const_iterator it;
                for (it = jacSparsity[i].begin(); it != jacSparsity[i].end(); ++it) {
                    size_t j = *it;
                    t[j] = true;
                }
            }
        }

        return true;
    }

    // ---------------------------------------------------------------------
    // Declare the AD<double> routine atomic_usead(id, ax, ay)
    CPPAD_USER_ATOMIC(atomic_use_lib,
                      CppAD::vector,
                      double,
                      atomic_use_lib_forward,
                      atomic_use_lib_reverse,
                      atomic_use_lib_for_jac_sparse,
                      atomic_use_lib_rev_jac_sparse,
                      atomic_use_lib_rev_hes_sparse)

}

using namespace CppAD;
using namespace std;

TEST_F(CppADCGDynamicAtomicTest, DynamicForRev) {
    using namespace std;
    using CppAD::vector;

    typedef double Base;
    typedef CppAD::CG<Base> CGD;
    typedef CppAD::AD<CGD> ADCG;

    vector< AD<double> > ax(n);
    for (size_t j = 0; j < n; j++)
        ax[j] = x[j];

    // declare independent variables and start tape recording
    CppAD::Independent(ax);

    vector<AD<double> > ay(m);

    // call user function and store atomic_use_lib(x) in au[0] 
    size_t id = 0; // not used
    CppAD::atomic_use_lib(id, ax, ay);

    // create f2: x -> y and stop tape recording
    ADFun<double> f2(ax, ay);

    /**
     * Test zero order
     */
    vector<CGD> xOrig(n);
    for (size_t j = 0; j < n; j++)
        xOrig[j] = x[j];

    vector<CGD> yOrig = _fun->Forward(0, xOrig);
    vector<double> yInner = CppADCGDynamicAtomicTest::MODEL->ForwardZero(x);
    vector<double> yOutter = f2.Forward(0, x);

    compareValues(yInner, yOrig);
    compareValues(yOutter, yOrig);

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
        vector<double> y_pInner = CppADCGDynamicAtomicTest::MODEL->ForwardOne(tx);
        vector<double> y_pOutter = f2.Forward(1, x_p);

        x_p[j] = 0;
        x_pOrig[j] = 0;
        tx[j * k1 + 1] = 0;

        compareValues(y_pInner, y_pOrig);
        compareValues(y_pOutter, y_pOrig);
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
        vector<double> dwInner = CppADCGDynamicAtomicTest::MODEL->ReverseOne(tx, ty, w);
        vector<double> dwOutter = f2.Reverse(1, w);

        w[i] = 0;
        wOrig[i] = 0;

        compareValues(dwInner, dwOrig);
        compareValues(dwOutter, dwOrig);
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
        vector<double> dwInner = CppADCGDynamicAtomicTest::MODEL->ReverseTwo(tx, ty, py);
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
    compareValues(jacOutter, jacOrig);

    /**
     * Jacobian sparsity
     */
    const std::vector<bool> jacSparsityOrig = jacobianSparsity < std::vector<bool>, CGD > (*_fun);
    const std::vector<bool> jacSparsityOutter = jacobianSparsity < std::vector<bool>, double > (f2);

    compareBoolValues(jacSparsityOrig, jacSparsityOutter);

    /**
     * Sparse jacobian
     */
    jacOrig = _fun->SparseJacobian(xOrig);
    jacOutter = f2.SparseJacobian(x);
    compareValues(jacOutter, jacOrig);

    /**
     * Hessian
     */
    for (size_t i = 0; i < m; i++) {
        w[i] = 1;
        wOrig[i] = 1;
    }
    vector<CGD> hessOrig = _fun->Hessian(xOrig, wOrig);
    vector<double> hessOutter = f2.Hessian(x, w);
    compareValues(hessOutter, hessOrig);

    /**
     * Sparse Hessian
     */
    hessOrig = _fun->SparseHessian(xOrig, wOrig);
    hessOutter = f2.SparseHessian(x, w);
    compareValues(hessOutter, hessOrig);


    // -----------------------------------------------------------------
    // Free all temporary work space associated with user_atomic objects. 
    // (If there are future calls to user atomic functions, they will 
    // create new temporary work space.)
    CppAD::user_atomic<double>::clear();
}