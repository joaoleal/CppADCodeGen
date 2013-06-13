#ifndef CPPAD_CG_UTIL_INCLUDED
#define CPPAD_CG_UTIL_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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

namespace CppAD {

    template<class VectorBool, class Base>
    inline VectorBool jacobianForwardSparsity(ADFun<Base>& fun) {
        size_t n = fun.Domain();

        VectorBool r(n * n);
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < n; k++)
                r[j * n + k] = false;
            r[j * n + j] = true;
        }
        return fun.ForSparseJac(n, r);

    }

    template<class VectorBool, class Base>
    inline VectorBool jacobianReverseSparsity(ADFun<Base>& fun) {
        size_t m = fun.Range();

        VectorBool s(m * m);
        for (size_t i = 0; i < m; i++) {
            for (size_t k = 0; k < m; k++)
                s[i * m + k] = false;
            s[i * m + i] = true;
        }
        return fun.RevSparseJac(m, s);
    }

    template<class VectorSet, class Base>
    inline VectorSet jacobianForwardSparsitySet(ADFun<Base>& fun) {
        size_t n = fun.Domain();

        VectorSet r(n);
        for (size_t i = 0; i < n; i++)
            r[i].insert(i);

        return fun.ForSparseJac(n, r);
    }

    template<class VectorSet, class Base>
    inline VectorSet jacobianReverseSparsitySet(ADFun<Base>& fun) {
        size_t m = fun.Range();

        VectorSet s_s(m);
        for (size_t i = 0; i < m; i++)
            s_s[i].insert(i);

        return fun.RevSparseJac(m, s_s);
    }

    /**
     * Determines the Jacobian sparsity for a model
     * 
     * @param fun The model
     * @return The Jacobian sparsity
     */
    template<class VectorBool, class Base>
    inline VectorBool jacobianSparsity(ADFun<Base>& fun) {
        size_t m = fun.Range();
        size_t n = fun.Domain();

        if (n <= m) {
            // use forward mode 
            return jacobianForwardSparsity<VectorBool, Base> (fun);
        } else {
            // use reverse mode 
            return jacobianReverseSparsity<VectorBool, Base> (fun);
        }
    }

    /**
     * Determines the Jacobian sparsity for a model
     * 
     * @param fun The model
     * @return The Jacobian sparsity
     */
    template<class VectorSet, class Base>
    inline VectorSet jacobianSparsitySet(ADFun<Base>& fun) {
        size_t m = fun.Range();
        size_t n = fun.Domain();

        if (n <= m) {
            // use forward mode 
            return jacobianForwardSparsitySet<VectorSet, Base> (fun);
        } else {
            // use reverse mode 
            return jacobianReverseSparsitySet<VectorSet, Base> (fun);
        }
    }

    /**
     * Determines the sum of the hessian sparsities for all the dependent 
     * variables in a model
     * 
     * @param fun The model
     * @return The sum of the hessian sparsities
     */
    template<class VectorBool, class Base>
    inline VectorBool hessianSparsity(ADFun<Base>& fun) {
        size_t m = fun.Range();
        size_t n = fun.Domain();

        /**
         * Determine the sparsity pattern p for Hessian of w^T F
         */
        std::vector<bool> r(n * n); // identity matrix
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < n; k++)
                r[j * n + k] = false;
            r[j * n + j] = true;
        }
        fun.ForSparseJac(n, r);

        std::vector<bool> s(m, true);
        return fun.RevSparseHes(n, s);
    }

    /**
     * Determines the hessian sparsity for a given dependent variable/equation
     * in a model
     * 
     * @param fun The model
     * @param i The dependent variable/equation index
     * @return The hessian sparsity
     */
    template<class VectorBool, class Base>
    inline VectorBool hessianSparsity(ADFun<Base>& fun, size_t i) {
        size_t m = fun.Range();
        size_t n = fun.Domain();

        /**
         * Determine the sparsity pattern p for Hessian of w^T F
         */
        std::vector<bool> r(n * n); // identity matrix
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < n; k++)
                r[j * n + k] = false;
            r[j * n + j] = true;
        }
        fun.ForSparseJac(n, r);

        std::vector<bool> s(m, false);
        s[i] = true;
        return fun.RevSparseHes(n, s);
    }

    template<class VectorBool>
    inline void generateSparsityIndexes(const VectorBool& sparsity,
                                        size_t m,
                                        size_t n,
                                        std::vector<size_t>& row,
                                        std::vector<size_t>& col) {
        assert(sparsity.size() == m * n);

        // determine total number of non zeros
        size_t nnz = 0;
        for (size_t i = 0; i < sparsity.size(); i++) {
            if (sparsity[i])
                nnz++;
        }

        row.resize(nnz);
        col.resize(nnz);

        // save the indexes
        nnz = 0;
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                if (sparsity[i * n + j]) {
                    row[nnz] = i;
                    col[nnz] = j;
                    nnz++;
                }
            }
        }

        assert(nnz == row.size());
    }

    inline void generateSparsityIndexes(const std::vector< std::set<size_t> >& sparsity,
                                        std::vector<size_t>& row,
                                        std::vector<size_t>& col) {
        std::vector< std::set<size_t> >::const_iterator rowIt;

        // determine total number of non zeros
        size_t nnz = 0;
        for (rowIt = sparsity.begin(); rowIt != sparsity.end(); ++rowIt) {
            nnz += rowIt->size();
        }

        row.resize(nnz);
        col.resize(nnz);

        // save the indexes
        nnz = 0;
        size_t i = 0;
        for (rowIt = sparsity.begin(); rowIt != sparsity.end(); ++rowIt, i++) {
            size_t rownnz = rowIt->size();
            std::fill(row.begin() + nnz, row.begin() + nnz + rownnz, i);
            std::copy(rowIt->begin(), rowIt->end(), col.begin() + nnz);
            nnz += rownnz;
        }
    }

    template<class VectorBool, class Base>
    void zeroOrderDependency(ADFun<Base>& fun,
                             const VectorBool& vx,
                             VectorBool& vy) {
        size_t m = fun.Range();
        CPPADCG_ASSERT_KNOWN(vx.size() >= fun.Domain(), "Invalid vx size");
        CPPADCG_ASSERT_KNOWN(vy.size() >= m, "Invalid vy size");
        
        typedef vector<std::set<size_t> > VectorSet;
        
        const VectorSet jacSparsity = jacobianSparsitySet<VectorSet, Base>(fun);
        
        std::set<size_t>::const_iterator it;
        for (size_t i = 0; i < m; i++) {
            for (it = jacSparsity[i].begin(); it != jacSparsity[i].end(); ++it) {
                size_t j = *it;
                if (vx[j]) {
                    vy[i] = true;
                    break;
                }
            }
        }
    }

    template<class VectorSet>
    inline bool isIdentityPattern(const VectorSet& pattern, size_t mRows) {
        assert(pattern.size() >= mRows);

        for (size_t i = 0; i < mRows; i++) {
            if (pattern[i].size() != 1 || *pattern[i].begin() != i) {
                return false;
            }
        }
        return true;
    }

    template<class VectorSet>
    inline VectorSet transposePattern(const VectorSet& pattern, size_t mRows, size_t nCols) {
        assert(pattern.size() >= mRows);

        VectorSet transpose(nCols);
        for (size_t i = 0; i < mRows; i++) {
            std::set<size_t>::const_iterator it;
            for (it = pattern[i].begin(); it != pattern[i].end(); ++it) {
                transpose[*it].insert(i);
            }
        }
        return transpose;
    }

    template<class VectorSet, class VectorSet2>
    inline void transposePattern(const VectorSet& pattern, size_t mRows, VectorSet2& transpose) {
        assert(pattern.size() >= mRows);

        for (size_t i = 0; i < mRows; i++) {
            std::set<size_t>::const_iterator it;
            for (it = pattern[i].begin(); it != pattern[i].end(); ++it) {
                transpose[*it].insert(i);
            }
        }
    }

    /**
     * Computes the resulting sparsity from adding one matrix to another:
     * R += A
     * 
     * @param a The matrix to be added to the result
     * @param result the resulting sparsity matrix
     */
    template<class VectorSet, class VectorSet2>
    inline void addMatrixSparsity(const VectorSet& a,
                                  VectorSet2& result) {
        assert(result.size() == a.size());

        for (size_t i = 0; i < a.size(); i++) {
            result[i].insert(a[i].begin(), a[i].end());
        }
    }

    /**
     * Computes the resulting sparsity from the multiplying of two matrices:
     * R += A * B
     * 
     * @param a The left matrix in the multiplication
     * @param b The right matrix in the multiplication
     * @param result the resulting sparsity matrix
     * @param m The number of rows of A
     * @param n The number of columns of A and rows of B
     * @param q The number of columns of B and the result
     */
    template<class VectorSet, class VectorSet2>
    inline void multMatrixMatrixSparsity(const VectorSet& a,
                                         const VectorSet2& b,
                                         CppAD::vector< std::set<size_t> >& result,
                                         size_t m,
                                         size_t n,
                                         size_t q) {
        assert(a.size() >= m);
        assert(b.size() >= n);
        assert(result.size() >= m);

        //check if b is identity
        if (n == q) {
            if (isIdentityPattern(b, n)) {
                for (size_t i = 0; i < m; i++) {
                    result[i] = a[i];
                }
                return;
            }
        }

        VectorSet2 bb = transposePattern(b, n, q);

        for (size_t jj = 0; jj < q; jj++) { //loop columns of b
            for (size_t i = 0; i < m; i++) {
                std::set<size_t>::const_iterator it;
                for (it = a[i].begin(); it != a[i].end(); ++it) {
                    if (bb[jj].find(*it) != bb[jj].end()) {
                        result[i].insert(jj);
                        break;
                    }
                }
            }
        }
    }

    /**
     * Computes the resulting sparsity from multiplying two matrices:
     * R += A^T * B
     * 
     * @param a The left matrix in the multiplication
     * @param b The right matrix in the multiplication
     * @param result the resulting sparsity matrix
     * @param m The number of rows of A and rows of B
     * @param n The number of columns of A and rows of the result
     * @param q The number of columns of B and the result
     */
    template<class VectorSet, class VectorSet2>
    inline void multMatrixTransMatrixSparsity(const VectorSet& a,
                                              const VectorSet2& b,
                                              CppAD::vector< std::set<size_t> >& result,
                                              size_t m,
                                              size_t n,
                                              size_t q) {
        assert(a.size() >= m);
        assert(b.size() >= m);
        assert(result.size() >= n);

        //check if B is empty
        bool empty = true;
        for (size_t i = 0; i < m; i++) {
            if (b[i].size() > 0) {
                empty = false;
                break;
            }
        }
        if (empty) {
            return; //nothing to do
        }

        //check if A is identity
        if (m == n && isIdentityPattern(a, m)) {
            for (size_t i = 0; i < n; i++) {
                result[i] = b[i];
            }
            return;
        }

        //check if B is identity
        if (m == q && isIdentityPattern(b, m)) {
            transposePattern(a, m, result);
            return;
        }

        VectorSet aa = transposePattern(a, m, n);
        VectorSet2 bb = transposePattern(b, m, q);

        for (size_t jj = 0; jj < q; jj++) { //loop columns of b
            for (size_t i = 0; i < n; i++) {
                std::set<size_t>::const_iterator it;
                for (it = aa[i].begin(); it != aa[i].end(); ++it) {
                    if (bb[jj].find(*it) != bb[jj].end()) {
                        result[i].insert(jj);
                        break;
                    }
                }
            }
        }
    }

    /**
     * Computes the transpose of the resulting sparsity from multiplying two
     * matrices:
     * (R += A * B)^T
     * 
     * @param a The TRANSPOSE of the left matrix in the multiplication
     * @param b The right matrix in the multiplication
     * @param result the TRANSPOSE of the resulting sparsity matrix
     * @param m The number of rows of B
     * @param n The number of columns of B
     * @param q The number of rows of A and the result
     */
    template<class VectorSet, class VectorSet2>
    inline void multMatrixMatrixSparsityTrans(const VectorSet& aT,
                                              const VectorSet2& b,
                                              CppAD::vector< std::set<size_t> >& rT,
                                              size_t m,
                                              size_t n,
                                              size_t q) {
        assert(aT.size() >= m);
        assert(b.size() >= m);

        //check if b is empty
        bool empty = true;
        for (size_t i = 0; i < m; i++) {
            if (b[i].size() > 0) {
                empty = false;
                break;
            }
        }
        if (empty) {
            return; //nothing to do
        }

        //check if a is identity
        if (m == q && isIdentityPattern(aT, m)) {
            transposePattern(b, m, rT);
            return;
        }

        VectorSet a = transposePattern(aT, m, q);
        VectorSet2 bT = transposePattern(b, m, n);

        for (size_t jj = 0; jj < n; jj++) { //loop columns of b
            for (size_t i = 0; i < q; i++) {
                std::set<size_t>::const_iterator it;
                for (it = a[i].begin(); it != a[i].end(); ++it) {
                    if (bT[jj].find(*it) != bT[jj].end()) {
                        rT[jj].insert(i);
                        break;
                    }
                }
            }
        }
    }
}

#endif

