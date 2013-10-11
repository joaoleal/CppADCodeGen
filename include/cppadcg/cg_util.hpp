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
    inline VectorBool hessianSparsity(ADFun<Base>& fun, bool transpose = false) {
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
        return fun.RevSparseHes(n, s, transpose);
    }

    template<class VectorSet, class Base>
    inline VectorSet hessianSparsitySet(ADFun<Base>& fun, bool transpose = false) {
        size_t m = fun.Range();

        std::set<size_t> w;
        for (size_t i = 0; i < m; i++) {
            w.insert(i);
        }
        return hessianSparsitySet<VectorSet, Base>(fun, w, transpose);
    }

    template<class VectorSet, class Base>
    inline VectorSet hessianSparsitySet(ADFun<Base>& fun, const std::set<size_t>& w, bool transpose = false) {
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

    /**
     * Determines the hessian sparsity for a given dependent variable/equation
     * in a model
     * 
     * @param fun The model
     * @param i The dependent variable/equation index
     * @return The hessian sparsity
     */
    template<class VectorBool, class Base>
    inline VectorBool hessianSparsity(ADFun<Base>& fun, size_t i, bool transpose = false) {
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
        return fun.RevSparseHes(n, s, transpose);
    }

    template<class VectorSet, class Base>
    inline VectorSet hessianSparsitySet(ADFun<Base>& fun, size_t i, bool transpose = false) {
        size_t n = fun.Domain();

        VectorSet r(n); // identity matrix
        for (size_t j = 0; j < n; j++)
            r[j].insert(j);
        fun.ForSparseJac(n, r);

        VectorSet s(1);
        s[0].insert(i);

        return fun.RevSparseHes(n, s, transpose);
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
    inline void transposePattern(const VectorSet& pattern, VectorSet2& transpose) {
        transposePattern(pattern, pattern.size(), transpose);
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
     * @param q The number of columns of B and the result
     */
    template<class VectorSet, class VectorSet2>
    inline void multMatrixMatrixSparsity(const VectorSet& a,
                                         const VectorSet2& b,
                                         CppAD::vector< std::set<size_t> >& result,
                                         size_t q) {
        multMatrixMatrixSparsity(a, b, result, a.size(), b.size(), q);
    }

    /**
     * Computes the resulting sparsity from the multiplying of two matrices:
     * R += A * B
     * 
     * Optimized for when the B matrix has less elements than A.
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

        VectorSet2 bt = transposePattern(b, n, q);
        std::set<size_t>::const_iterator itRowb;

        for (size_t jj = 0; jj < q; jj++) { //loop columns of b
            const std::set<size_t>& colB = bt[jj];
            if (colB.size() > 0) {
                for (size_t i = 0; i < m; i++) {
                    const std::set<size_t>& rowA = a[i];
                    for (itRowb = colB.begin(); itRowb != colB.end(); ++itRowb) {
                        if (rowA.find(*itRowb) != rowA.end()) {
                            result[i].insert(jj);
                            break;
                        }
                    }
                }
            }
        }
    }

    /**
     * Computes the resulting sparsity from multiplying two matrices:
     * R += A^T * B
     * 
     * Optimized for when the B matrix has less elements than A.
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

        VectorSet at = transposePattern(a, m, n);
        VectorSet2 bt = transposePattern(b, m, q);
        std::set<size_t>::const_iterator itRowb;

        for (size_t jj = 0; jj < q; jj++) { //loop columns of b
            const std::set<size_t>& colB = bt[jj];
            if (colB.size() > 0) {
                for (size_t i = 0; i < n; i++) {
                    const std::set<size_t>& rowAt = at[i];
                    if (rowAt.size() > 0) {
                        for (itRowb = colB.begin(); itRowb != colB.end(); ++itRowb) {
                            if (rowAt.find(*itRowb) != rowAt.end()) {
                                result[i].insert(jj);
                                break;
                            }
                        }
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

    template<class VectorBool>
    void printSparsityPattern(const VectorBool& sparsity,
                              const std::string& name,
                              size_t m, size_t n) {
        size_t width = std::ceil(std::log10((m > n) ? m : n));
        if (!name.empty()) {
            std::cout << name << "  sparsity:\n";
        }
        for (size_t i = 0; i < m; i++) {
            std::cout << " " << std::setw(width) << i << ": ";
            for (size_t j = 0; j < n; j++) {
                if (sparsity[i * n + j]) {
                    std::cout << std::setw(width) << j << " ";
                } else {
                    std::cout << std::setw(width) << " " << " ";
                }
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    template<class VectorSet>
    void printSparsityPattern(const VectorSet& sparsity,
                              const std::string& name) {
        size_t maxDim = sparsity.size();
        for (size_t i = 0; i < sparsity.size(); i++) {
            if (sparsity[i].size() > 0 && *sparsity[i].rbegin() > maxDim) {
                maxDim = *sparsity[i].rbegin();
            }
        }

        size_t width = std::ceil(std::log10(maxDim));
        if (!name.empty()) {
            std::cout << name << "  sparsity:\n";
        }
        std::set<size_t>::const_iterator it;
        for (size_t i = 0; i < sparsity.size(); i++) {
            std::cout << " " << std::setw(width) << i << ": ";
            long last = -1;
            for (it = sparsity[i].begin(); it != sparsity[i].end(); ++it) {
                size_t j = *it;
                if (j != 0 && long(j) != last + 1) {
                    std::cout << std::setw((j - last - 1) * (width + 1)) << " ";
                }
                std::cout << std::setw(width) << j << " ";
                last = j;
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    template<class VectorSize>
    void printSparsityPattern(const VectorSize& row,
                              const VectorSize& col,
                              const std::string& name,
                              size_t m) {
        vector<std::set<size_t> > sparsity(m);
        generateSparsitySet(row, col, sparsity);
        printSparsityPattern(sparsity, name);
    }

    inline bool intersects(const std::set<size_t>& a,
                           const std::set<size_t>& b) {
        if (a.empty() || b.empty()) {
            return false;
        } else if (*a.rbegin() < *b.begin() ||
                *a.begin() > *b.rbegin()) {
            return false;
        }

        if (a.size() < b.size()) {
            std::set<size_t>::const_iterator ita;
            for (ita = a.begin(); ita != a.end(); ++ita) {
                if (b.find(*ita) != b.end()) {
                    return true;
                }
            }
        } else {
            std::set<size_t>::const_iterator itb;
            for (itb = b.begin(); itb != b.end(); ++itb) {
                if (a.find(*itb) != a.end()) {
                    return true;
                }
            }
        }

        return false;
    }

    /**
     * Finds the first non-null code handler
     * 
     * @param ty The array to search in
     * @return The first code handler found or NULL if none was found
     */
    template<class Base>
    inline CodeHandler<Base>* findHandler(const vector<CG<Base> >& ty) {
        for (size_t i = 0; i < ty.size(); i++) {
            if (ty[i].getCodeHandler() != NULL) {
                return ty[i].getCodeHandler();
            }
        }
        return NULL;
    }

    template<class Base>
    inline Argument<Base> asArgument(const CG<Base>& tx) {
        if (tx.isParameter()) {
            return Argument<Base>(tx.getValue());
        } else {
            return Argument<Base>(*tx.getOperationNode());
        }
    }

    template<class Base>
    inline std::vector<Argument<Base> > asArguments(const vector<CG<Base> >& tx) {
        std::vector<Argument<Base> > arguments(tx.size());
        for (size_t i = 0; i < arguments.size(); i++) {
            arguments[i] = asArgument(tx[i]);
        }
        return arguments;
    }

    /**
     * Smart vector of pointers.
     * Deletes all vector values on destruction.
     */
    template<class Base>
    class SmartVectorPointer {
    public:
        std::vector<Base*> v;

        inline SmartVectorPointer() {
        }

        inline SmartVectorPointer(size_t size) :
            v(size) {
        }

        ~SmartVectorPointer() {
            for (size_t i = 0; i < v.size(); i++) {
                delete v[i];
            }
        }
    };

    /**
     * Smart set of pointers.
     * Deletes all set values on destruction.
     */
    template<class Base>
    class SmartSetPointer {
    public:
        std::set<Base*> s;

        inline SmartSetPointer() {
        }

        inline SmartSetPointer(const std::set<Base*>& s_) {
            s.swap(s_);
        }

        ~SmartSetPointer() {
            typename std::set<Base*>::const_iterator it;
            for (it = s.begin(); it != s.end(); ++it) {
                delete *it;
            }
        }
    };

    template<class Key, class Value>
    class SmartMapValuePointer {
    public:
        std::map<Key, Value*> m;

        std::map<Key, Value*> release() {
            std::map<Key, Value*> result;
            result.swap(m);
            return result;
        }

        ~SmartMapValuePointer() {
            typename std::map<Key, Value*>::const_iterator it;
            for (it = m.begin(); it != m.end(); ++it) {
                delete it->second;
            }
        }
    };

    /***************************************************************************
     * map related
     **************************************************************************/

    /**
     * Gets all the keys present in a map
     * 
     * @param map the map from which to get the keys from
     * @param keys the map keys will be inserted into this set
     */
    template<class Key, class Value>
    void mapKeys(const std::map<Key, Value>& map, std::set<Key>& keys) {
        typename std::map<Key, Value>::const_iterator it;
        for (it = map.begin(); it != map.end(); ++it) {
            keys.insert(keys.end(), it->first);
        }
    }

    /**
     * Gets all the keys present in a map
     * 
     * @param map the map from which to get the keys from
     * @param keys the map keys will be saved in this vector
     */
    template<class Key, class Value>
    void mapKeys(const std::map<Key, Value>& map, std::vector<Key>& keys) {
        keys.resize(map.size());

        size_t i = 0;
        typename std::map<Key, Value>::const_iterator it;
        for (it = map.begin(); it != map.end(); ++it, i++) {
            keys[i] = it->first;
        }
    }

    /**
     * Checks if a map has only a set of keys.
     * 
     * @param map The map 
     * @param keys The keys
     * @return true if all the keys and only these keys where found in the map
     */
    template<class Key, class Value>
    bool compareMapKeys(const std::map<Key, Value>& map, const std::set<Key>& keys) {
        if (map.size() != keys.size())
            return false;

        typename std::map<Key, Value>::const_iterator itm = map.begin();
        typename std::set<Key>::const_iterator itk = keys.begin();
        for (; itm != map.end(); ++itm, ++itk) {
            if (itm->first != *itk)
                return false;
        }

        return true;
    }

    /**
     * Creates a new map with only a given set of keys
     * 
     * @param m The map to be filtered
     * @param keys the keys (the filter) to be retrieved from the map 
     * @return a new map only with the keys found in provided filter
     */
    template<class Key, class Value>
    inline std::map<Key, Value> filterBykeys(const std::map<Key, Value>& m,
                                             const std::set<Key>& keys) {
        std::map<Key, Value> filtered;

        typename std::map<Key, Value>::const_iterator itM;

        typename std::set<Key>::const_iterator itK;
        for (itK = keys.begin(); itK != keys.end(); ++itK) {
            itM = m.find(*itK);
            if (itM != m.end()) {
                filtered[itM->first] = itM->second;
            }
        }
        return filtered;
    }

    /**
     * Compares two sets
     * 
     * @param s1 the first set
     * @param s2 the second set
     * @return -1 if the first set is considered lower than the second,
     *         0 if they have all the same elements
     *         1 if the second set is considered lower than the first.
     */
    template<class T>
    inline int compare(const std::set<T>& s1, const std::set<T>& s2) {
        if (s1.size() < s2.size()) {
            return -1;
        } else if (s1.size() > s2.size()) {
            return 1;
        } else {
            typename std::set<T>::const_iterator it1, it2;
            for (it1 = s1.begin(), it2 = s2.begin(); it1 != s1.end(); ++it1, ++it2) {
                if (*it1 < *it2) {
                    return -1;
                } else if (*it1 > *it2) {
                    return 1;
                }
            }
            return 0;
        }
    }

    template<class Base>
    inline void printModel(ADFun<CG<Base> >& fun) {
        std::vector<std::string> depNames;
        std::vector<std::string> indepNames;
        printModel(fun, depNames, indepNames);
    }

    template<class Base>
    inline void printModel(ADFun<CG<Base> >& fun,
                           const std::vector<std::string>& depNames,
                           const std::vector<std::string>& indepNames) {
        assert(depNames.size() <= fun.Range());
        assert(indepNames.size() <= fun.Domain());

        CodeHandler<Base> handler;

        vector<CG<Base> > indep0(fun.Domain());
        handler.makeVariables(indep0);

        vector<CG<Base> > dep0 = fun.Forward(0, indep0);

        CLanguage<double> langC("double");

        /**
         * generate the source code
         */
        CLangCustomVariableNameGenerator<double> nameGen(depNames, indepNames,
                                                         "y", "x", "z", "array");

        std::ostringstream code;
        handler.generateCode(code, langC, dep0, nameGen);
        std::cout << "\n" << code.str() << std::endl;
    }
}

#endif

