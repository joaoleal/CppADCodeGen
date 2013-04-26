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
            return jacobianForwardSparsity<VectorBool, Base > (fun);
        } else {
            // use reverse mode 
            return jacobianReverseSparsity<VectorBool, Base >(fun);
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
     * @param q The number of columns of b and the result
     */
    template<class VectorSet, class VectorSet2>
    inline void multMatrixMatrixSparsity(const VectorSet& a,
                                         const VectorSet2& b,
                                         CppAD::vector< std::set<size_t> >& result,
                                         size_t q) {
        const size_t m = a.size();
        const size_t n = b.size();
        assert(result.size() == m);

        //check if b is identity
        if (n == q) {
            bool identity = true;
            for (size_t i = 0; i < n; i++) {
                if (b[i].size() != 1 || *b[i].begin() != i) {
                    identity = false;
                    break;
                }
            }
            if (identity) {
                for (size_t i = 0; i < m; i++) {
                    result[i] = a[i];
                }
                return;
            }
        }

        VectorSet2 bb(q);
        for (size_t i = 0; i < n; i++) {
            std::set<size_t>::const_iterator it;
            for (it = b[i].begin(); it != b[i].end(); ++it) {
                bb[*it].insert(i);
            }
        }

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
     * @param q The number of columns of b and the result
     */
    template<class VectorSet, class VectorSet2>
    inline void multMatrixTransMatrixSparsity(const VectorSet& a,
                                              const VectorSet2& b,
                                              CppAD::vector< std::set<size_t> >& result,
                                              size_t q) {
        const size_t m = a.size();
        const size_t n = b.size();
        assert(m == n);

        //check if b is empty
        bool empty = true;
        for (size_t i = 0; i < n; i++) {
            if (b[i].size() > 0) {
                empty = false;
                break;
            }
        }
        if (empty) {
            return; //nothing to do
        }

        //check if b is identity
        if (n == q) {
            bool identity = true;
            for (size_t i = 0; i < n; i++) {
                if (b[i].size() != 1 || *b[i].begin() != i) {
                    identity = false;
                    break;
                }
            }
            if (identity) {
                for (size_t i = 0; i < m; i++) {
                    std::set<size_t>::const_iterator it;
                    for (it = a[i].begin(); it != a[i].end(); ++it) {
                        result[*it].insert(i);
                    }
                }
                return;
            }
        }

        VectorSet2 bb(q);
        for (size_t i = 0; i < n; i++) {
            std::set<size_t>::const_iterator it;
            for (it = b[i].begin(); it != b[i].end(); ++it) {
                bb[*it].insert(i);
            }
        }

        VectorSet2 aa(n);
        for (size_t i = 0; i < m; i++) {
            std::set<size_t>::const_iterator it;
            for (it = a[i].begin(); it != a[i].end(); ++it) {
                aa[*it].insert(i);
            }
        }

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
}

#endif

