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

#include <vector>
#include <valarray>

#include <cppadcg/cg.hpp>
#include <gtest/gtest.h>
#include "CppADCGTest.hpp"

using namespace CppAD;
using namespace std;

class SparseJacHes : public CppADCGTest {
public:
    typedef std::vector< std::set<size_t> > std_vector_set;
    typedef CppAD::vector< std::set<size_t> > cppad_vector_set;
public:

    template <class VectorBase, class VectorSet, class VectorSize>
    void set_case() {

        size_t i, j;

        // domain space vector
        size_t n = 3;
        CPPAD_TESTVECTOR(AD<double>) X(n);
        for (i = 0; i < n; i++)
            X[i] = AD<double> (0);

        // declare independent variables and starting recording
        Independent(X);

        size_t m = 2;
        CPPAD_TESTVECTOR(AD<double>) Y(m);
        Y[0] = X[0] * X[1] * X[2] + X[0];
        Y[1] = X[1] * X[2] + 10;

        // create f: x -> y and stop tape recording
        ADFun<double> f(X, Y);

        // new value for the independent variable vector
        VectorBase x(n);
        for (i = 0; i < n; i++)
            x[i] = 3 * double(i + 1);

        // second derivative of y[1] 
        VectorBase w(m);
        w[0] = 1.;
        w[1] = 2.;

        /**
         * zero-order
         */
        VectorBase check_y = f.Forward(0, x);

        /**
         * Jacobian
         */
        // sparsity pattern
        VectorSet s(m), p(m);
        for (i = 0; i < m; i++)
            s[i].insert(i);
        VectorSet jsparsity = f.RevSparseJac(m, s);

        // evaluate
        size_t jnnz = 2;
        VectorSize jrow(jnnz), jcol(jnnz);
        jrow[0] = 0;
        jcol[0] = 2;
        jrow[1] = 1;
        jcol[1] = 1;
        VectorBase check_jac(jnnz);
        sparse_jacobian_work jwork;

        size_t n_sweep_j = f.SparseJacobianForward(x, jsparsity, jrow, jcol, check_jac, jwork);

        /**
         * Hessian
         */
        // determine the sparsity pattern p for Hessian of w^T F
        VectorSet r(n);
        for (j = 0; j < n; j++)
            r[j].insert(j);
        f.ForSparseJac(n, r);
        //
        s.resize(1);
        for (i = 0; i < m; i++)
            if (w[i] != 0)
                s[0].insert(i);
        VectorSet hsparsity = f.RevSparseHes(n, s);

        // evaluate
        size_t hnnz = 5;
        VectorSize hrow(hnnz), hcol(hnnz);
        hrow[0] = 0;
        hcol[0] = 1;
        hrow[1] = 0;
        hcol[1] = 2;
        hrow[2] = 1;
        hcol[2] = 0;
        hrow[3] = 1;
        hcol[3] = 2;
        hrow[4] = 2;
        hcol[4] = 1;
        VectorBase check_hes(hnnz);
        sparse_hessian_work hwork;

        size_t n_sweep_h = f.SparseHessian(x, w, hsparsity, hrow, hcol, check_hes, hwork);

        size_t n_sweep = 1 + n_sweep_j + 2 * n_sweep_h;

        /**
         * Jacobian + Hessian
         */
        VectorBase y(m);
        VectorBase jac(jnnz);
        VectorBase hes(hnnz);

        SparseForjacHessianWork work;
        size_t n_sweep_jh = 1 + SparseForJacHessian(f, x, w,
                                                    y,
                                                    jsparsity, jrow, jcol, jac,
                                                    hsparsity, hrow, hcol, hes,
                                                    work);
        ASSERT_TRUE(compareValues(y, check_y));
        ASSERT_TRUE(compareValues(jac, check_jac));
        ASSERT_TRUE(compareValues(hes, check_hes));

    }

};

TEST_F(SparseJacHes, CppADVectorVectorSet) {
    this->set_case< CppAD::vector<double>, std_vector_set, std::vector<size_t> >();
}

TEST_F(SparseJacHes, valarrayVectorSet) {
    set_case< std::valarray<double>, std_vector_set, std::vector<size_t> >();
}

TEST_F(SparseJacHes, VectorCppADVectorSet) {
    set_case< std::vector<double>, cppad_vector_set, std::vector<size_t> >();
}

TEST_F(SparseJacHes, CppADVectorCppADVectorSet) {
    set_case< CppAD::vector<double>, cppad_vector_set, std::vector<size_t> >();
}
