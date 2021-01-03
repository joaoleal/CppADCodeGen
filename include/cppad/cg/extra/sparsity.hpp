#ifndef CPPAD_CG_SPARSITY_INCLUDED
#define CPPAD_CG_SPARSITY_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *    Copyright (C) 2019 Joao Leal
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

namespace CppAD {
namespace cg {

template<class Base>
inline sparse_rc<CppAD::vector<size_t>> jacobianForwardSparsity(ADFun<Base>& fun,
                                                                bool internal_bool=true) {
    bool transpose = false;
    bool dependency = true;
    size_t n = fun.Domain();

    sparse_rc<CppAD::vector<size_t> > pattern_in;
    pattern_in.resize(n, n, n);

    for (size_t k = 0; k < n; ++k)
        pattern_in.set(k, k, k);

    sparse_rc<CppAD::vector<size_t> > jac_sparsity;

    fun.for_jac_sparsity(pattern_in, transpose, dependency, internal_bool, jac_sparsity);

    if (internal_bool)
        fun.size_forward_bool(0);
    else
        fun.size_forward_set(0);

    return jac_sparsity;
}

template<class Base>
inline sparse_rc<CppAD::vector<size_t>> jacobianReverseSparsity(ADFun<Base>& fun,
                                                                bool internal_bool=true) {
    bool transpose = false;
    bool dependency = true;
    size_t m = fun.Range();

    sparse_rc<CppAD::vector<size_t> > pattern_in;
    pattern_in.resize(m, m, m);

    for (size_t k = 0; k < m; ++k)
        pattern_in.set(k, k, k);

    sparse_rc<CppAD::vector<size_t> > jac_sparsity;

    fun.rev_jac_sparsity(pattern_in, transpose, dependency, internal_bool, jac_sparsity);

    if (internal_bool)
        fun.size_forward_bool(0);
    else
        fun.size_forward_set(0);

    return jac_sparsity;
}

template<class VectorBool, class Base>
inline VectorBool jacobianForwardSparsityBool(ADFun<Base>& fun) {
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
inline VectorBool jacobianReverseSparsityBool(ADFun<Base>& fun) {
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
inline VectorBool jacobianSparsityBool(ADFun<Base>& fun) {
    size_t m = fun.Range();
    size_t n = fun.Domain();

    if (n <= m) {
        // use forward mode 
        return jacobianForwardSparsityBool<VectorBool, Base> (fun);
    } else {
        // use reverse mode 
        return jacobianReverseSparsityBool<VectorBool, Base> (fun);
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

template<class Base>
inline sparse_rc<CppAD::vector<size_t>> jacobianSparsity(ADFun<Base>& fun,
                                                         bool internal_bool=true) {
    size_t n = fun.Domain();
    size_t m = fun.Range();

    if (n <= m) {
        // use forward mode
        return jacobianForwardSparsity(fun, internal_bool);
    } else {
        // use reverse mode
        return jacobianReverseSparsity(fun, internal_bool);
    }
}

/**
 * Estimates the work load of forward vs reverse mode for the evaluation of
 * a Jacobian
 * 
 * @return true if the foward mode should be used, false for the reverse mode
 */
inline bool estimateBestJacobianADMode(const std::vector<size_t>& jacRows,
                                       const std::vector<size_t>& jacCols) {
    std::set<size_t> rows, cols;
    rows.insert(jacRows.begin(), jacRows.end());
    size_t workReverse = rows.size();
    cols.insert(jacCols.begin(), jacCols.end());
    size_t workForward = cols.size();

    return workForward <= workReverse;
}

/**
 * Determines the sum of the hessian sparsities for all the dependent 
 * variables in a model
 * 
 * @param fun The model
 * @return The sum of the hessian sparsities
 */
template<class VectorBool, class Base>
inline VectorBool hessianSparsityBool(ADFun<Base>& fun,
                                      bool transpose = false) {
    size_t m = fun.Range();
    size_t n = fun.Domain();

    /**
     * Determine the sparsity pattern p for Hessian of w^T F
     */
    VectorBool r(n * n); // identity matrix
    for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < n; k++)
            r[j * n + k] = false;
        r[j * n + j] = true;
    }
    fun.ForSparseJac(n, r);

    VectorBool s(m);
    for (size_t i = 0; i < m; i++)
        s[i] = true;
    return fun.RevSparseHes(n, s, transpose);
}

template<class VectorSet, class Base>
inline VectorSet hessianSparsitySet(ADFun<Base>& fun,
                                    const std::set<size_t>& w,
                                    bool transpose = false) {
    size_t n = fun.Domain();

    /**
     * Determine the sparsity pattern p for Hessian of w^T F
     */
    VectorSet r(n); // identity matrix
    for (size_t j = 0; j < n; j++)
        r[j].insert(j);
    fun.ForSparseJac(n, r);

    VectorSet s(1);
    s[0] = w;

    return fun.RevSparseHes(n, s, transpose);
}

template<class VectorSet, class Base>
inline VectorSet hessianSparsitySet(ADFun<Base>& fun,
                                    bool transpose = false) {
    size_t m = fun.Range();

    std::set<size_t> w;
    for (size_t i = 0; i < m; i++) {
        w.insert(i);
    }
    return hessianSparsitySet<VectorSet, Base>(fun, w, transpose);
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
inline VectorBool hessianSparsityBool(ADFun<Base>& fun,
                                      size_t i,
                                      bool transpose = false) {
    size_t m = fun.Range();
    size_t n = fun.Domain();

    /**
     * Determine the sparsity pattern p for Hessian of w^T F
     */
    VectorBool r(n * n); // identity matrix
    for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < n; k++)
            r[j * n + k] = false;
        r[j * n + j] = true;
    }
    fun.ForSparseJac(n, r);

    VectorBool s(m);
    for (size_t ii = 0; ii < m; ii++)
        s[ii] = false;
    s[i] = true;
    return fun.RevSparseHes(n, s, transpose);
}

template<class VectorSet, class Base>
inline VectorSet hessianSparsitySet(ADFun<Base>& fun,
                                    size_t i,
                                    bool transpose = false) {
    size_t n = fun.Domain();

    VectorSet r(n); // identity matrix
    for (size_t j = 0; j < n; j++)
        r[j].insert(j);
    fun.ForSparseJac(n, r);

    VectorSet s(1);
    s[0].insert(i);

    return fun.RevSparseHes(n, s, transpose);
}

template<class Base>
inline sparse_rc<CppAD::vector<size_t>> hessianSparsity(ADFun<Base>& fun,
                                                        bool internal_bool=true) {
    size_t m = fun.Range();
    size_t n = fun.Domain();

    CppAD::vector<bool> select_y(m), select_x(n);

    sparse_rc<CppAD::vector<size_t> > hes_sparsity;

    for (size_t i = 0; i < m; ++i)
        select_y[i] = true;

    if (n <= m) {
        for (size_t j = 0; j < n; ++j)
            select_x[j] = true;
        fun.for_hes_sparsity(select_x, select_y, internal_bool, hes_sparsity);
    } else {

        // sparsity pattern for the identity matrix
        sparse_rc<CppAD::vector<size_t> > pattern_in(n, n, n);
        for (size_t k = 0; k < n; k++) {
            pattern_in.set(k, k, k);
        }

        // forward Jacobian sparsity is stored in fun
        bool transpose = false;
        bool dependency = false;
        sparse_rc<CppAD::vector<size_t> > jac_pattern;
        fun.for_jac_sparsity(pattern_in, transpose, dependency, internal_bool, jac_pattern);

        fun.rev_hes_sparsity(select_y, transpose, internal_bool, hes_sparsity);
    }

    if (internal_bool)
        fun.size_forward_bool(0);
    else
        fun.size_forward_set(0);

    return hes_sparsity;
}

template<class VectorBool, class VectorSize>
inline void generateSparsityIndexes(const VectorBool& sparsity,
                                    size_t m,
                                    size_t n,
                                    VectorSize& row,
                                    VectorSize& col) {
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

template<class VectorSet, class VectorSize>
inline void generateSparsityIndexes(const VectorSet& sparsity,
                                    VectorSize& row,
                                    VectorSize& col) {
    size_t m = sparsity.size();

    // determine total number of non zeros
    size_t nnz = 0;
    for (size_t i = 0; i < m; i++) {
        nnz += sparsity[i].size();
    }

    row.resize(nnz);
    col.resize(nnz);
    if (nnz == 0)
        return;

    // save the indexes
    nnz = 0;
    for (size_t i = 0; i < m; i++) {
        const std::set<size_t>& rowSparsity = sparsity[i];
        size_t rowNnz = rowSparsity.size();
        std::fill(&row[0] + nnz, &row[0] + nnz + rowNnz, i);
        std::copy(rowSparsity.begin(), rowSparsity.end(), &col[0] + nnz);
        nnz += rowNnz;
    }
}

template<class VectorSet, class VectorSize>
inline void generateSparsitySet(const VectorSize& row,
                                const VectorSize& col,
                                VectorSet& sparsity) {
    assert(row.size() == col.size());

    size_t nnz = row.size();
    for (size_t e = 0; e < nnz; e++) {
        sparsity[row[e]].insert(col[e]);
    }
}

template<class VectorSet, class VectorSize>
inline void generateSparsitySet(const sparse_rc<VectorSize>& from,
                                VectorSet& to) {

    to.resize(from.nr());
    for (size_t k = 0; k < to.size(); ++k) {
        to[k].clear();
    }

    size_t nnz = from.nnz();
    const auto& row = from.row();
    const auto& col = from.col();

    for (size_t e = 0; e < nnz; e++) {
        to[row[e]].insert(col[e]);
    }
}

template<class VectorSet, class VectorSize>
inline void generateSparsity(const VectorSet& from,
                             size_t n,
                             sparse_rc<VectorSize>& to) {

    CPPADCG_ASSERT_KNOWN(to.nr() == from.size(), "Invalid sizes")

    size_t m = from.size();

    // determine total of non-zeros
    size_t nnz = 0;
    for (size_t i = 0; i < m; ++i)
        nnz += from[i].size();

    // reserve space
    to.resize(m, n, nnz);

    // copy indices
    size_t k = 0;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j: from[i]) {
            to.set(k, i, j);
            k++;
        }
    }
}

inline bool isAllTrue(ArrayView<const bool> select_x) noexcept {
    for (bool i : select_x)
        if (!i)
            return false;
    return true;
}

/**
 * See atomic_three<Base>::for_type
 */
inline void for_type(const sparse_rc<CppAD::vector<size_t>>& jac_pattern,
                     ArrayView<const ad_type_enum> type_x,
                     ArrayView<ad_type_enum> type_y) {

    size_t nr = jac_pattern.nr();
    size_t nnz = jac_pattern.nnz();
    const auto& row = jac_pattern.row();
    const auto& col = jac_pattern.col();

    CPPADCG_ASSERT_KNOWN(jac_pattern.nr() == type_y.size(), "Invalid type_y size")
    CPPADCG_ASSERT_KNOWN(jac_pattern.nc() == type_x.size(), "Invalid type_x size")

    // initialize type_y as constant_enum
    for (size_t i = 0; i < nr; ++i)
        type_y[i] = constant_enum;

    // loop over entries in Dependency pattern
    for (size_t k = 0; k < nnz; ++k) {
        size_t i = row[k];
        size_t j = col[k];
        type_y[i] = std::max(type_y[i], type_x[j]);
    }
}

/**
 * See atomic_three<Base>::for_type
 */
template<class Base>
void for_type(ADFun<Base>& fun,
              ArrayView<const ad_type_enum> type_x,
              ArrayView<ad_type_enum> type_y) {
    const sparse_rc<CppAD::vector<size_t>>& jac_pattern = CppAD::cg::jacobianSparsity(fun);

    for_type(jac_pattern, type_x, type_y);
}

/**
 * See atomic_three<Base>::rev_depend
 */
inline void rev_depend(const sparse_rc<CppAD::vector<size_t>>& jac_pattern,
                       ArrayView<const ad_type_enum> type_x,
                       ArrayView<bool> depend_x,
                       ArrayView<const bool> depend_y) {

    size_t nc = jac_pattern.nc();
    size_t nnz = jac_pattern.nnz();
    const auto& row = jac_pattern.row();
    const auto& col = jac_pattern.col();

    CPPAD_ASSERT_UNKNOWN(jac_pattern.nr() == depend_y.size())
    CPPAD_ASSERT_UNKNOWN(jac_pattern.nc() == depend_x.size())

    // initialize depend_x as false
    for (size_t j = 0; j < nc; ++j)
        depend_x[j] = false;

    // loop over entries in Dependency pattern
    for (size_t k = 0; k < nnz; ++k) {
        size_t i = row[k];
        size_t j = col[k];
        if (depend_y[i])
            depend_x[j] = true;
    }
}

template<class Base>
void rev_depend(ADFun<Base>& fun,
                ArrayView<const ad_type_enum> type_x,
                ArrayView<bool> depend_x,
                ArrayView<const bool> depend_y) {
    const sparse_rc<CppAD::vector<size_t>>& jac_pattern = CppAD::cg::jacobianSparsity(fun);
    rev_depend(jac_pattern, type_x, depend_x, depend_y);
}

inline void filter(const sparse_rc<CppAD::vector<size_t>>& from,
                   ArrayView<const bool> select_y,
                   ArrayView<const bool> select_x,
                   sparse_rc<CppAD::vector<size_t>>& to) {

    if (isAllTrue(select_x) && isAllTrue(select_y)) {
        to = from;
        return;
    }

    size_t fullNnz = from.nnz();
    auto& row = from.row();
    auto& col = from.col();

    //sparse_rc<CppAD::vector<size_t>> filtered(sparsity.nr(), sparsity.nc(), fullNnz);
    to.resize(from.nr(), from.nc(), fullNnz);

    size_t nnz = 0;

    for (size_t k = 0; k < fullNnz; ++k) {
        size_t i = row[k];
        size_t j = col[k];
        if (select_y[i] && select_x[j]) {
            to.set(nnz, i, j);
            nnz++;
        }
    }

    if (fullNnz != nnz) {
        to.resize(to.nr(), to.nc(), nnz);
    }
}

inline void filter(sparse_rc<CppAD::vector<size_t>>& sparsity,
                   ArrayView<const bool> select_y,
                   ArrayView<const bool> select_x) {

    sparse_rc<CppAD::vector<size_t>> filtered;
    filter(sparsity, select_y, select_x, filtered);


    if (sparsity.nnz() != filtered.nnz()) {
        sparsity = filtered;
    }
}

inline void filter(const sparse_rc<CppAD::vector<size_t>>& from,
                   ArrayView<const bool> select_x,
                   sparse_rc<CppAD::vector<size_t>>& to) {

    CPPADCG_ASSERT_UNKNOWN(from.nr() == from.nc())

    if (isAllTrue(select_x)) {
        to = from;
        return;
    }

    size_t fullNnz = from.nnz();
    auto& row = from.row();
    auto& col = from.col();

    to.resize(from.nr(), from.nc(), fullNnz);

    size_t nnz = 0;
    for (size_t k = 0; k < fullNnz; ++k) {
        size_t i = row[k];
        size_t j = col[k];
        if (select_x[i] && select_x[j]) {
            to.set(nnz, i, j);
            nnz++;
        }
    }

    if (fullNnz != nnz) {
        to.resize(to.nr(), to.nc(), nnz);
    }
}

inline void filter(sparse_rc<CppAD::vector<size_t>>& sparsity,
                   ArrayView<const bool> select_x) {

    sparse_rc<CppAD::vector<size_t>> filtered;
    filter(sparsity, select_x, filtered);

    if (sparsity.nnz() != filtered.nnz()) {
        sparsity = filtered;
    }
}

inline void filter(const std::vector<std::set<size_t> >& from,
                   ArrayView<const bool> filter,
                   sparse_rc<CppAD::vector<size_t>>& to) {

    CPPADCG_ASSERT_UNKNOWN(to.nr() == to.nc())
    CPPADCG_ASSERT_UNKNOWN(to.nr() == 0 || to.nr() == filter.size())
    CPPADCG_ASSERT_UNKNOWN(from.size() == filter.size())
    CPPADCG_ASSERT_UNKNOWN(to.nr() == from.size())

    if (isAllTrue(filter)) {
        generateSparsity(from, filter.size(), to);

    } else {
        size_t m = from.size();

        size_t fullNnz = 0;
        for (size_t i = 0; i < m; ++i)
            fullNnz += from[i].size();

        // reserve space
        to.resize(m, m, fullNnz);

        // copy indices
        size_t nnz = 0;
        for (size_t i = 0; i < m; ++i) {
            if (filter[i]) {
                for (size_t j: from[i]) {
                    if (filter[j]) {
                        to.set(nnz, i, j);
                        nnz++;
                    }
                }
            }
        }

        to.resize(m, m, nnz);
    }
}

} // END cg namespace
} // END CppAD namespace

#endif

