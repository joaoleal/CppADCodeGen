#ifndef CPPAD_CG_UTIL_INCLUDED
#define	CPPAD_CG_UTIL_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

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

    template<class Base>
    inline std::vector< std::set<size_t> > jacobianReverseSparsitySet(ADFun<Base>& fun) {
        size_t m = fun.Range();

        std::vector< std::set<size_t> > s_s(m);
        for (size_t i = 0; i < m; i++)
            s_s[i].insert(i);

        return fun.RevSparseJac(m, s_s);
    }

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

}

#endif

