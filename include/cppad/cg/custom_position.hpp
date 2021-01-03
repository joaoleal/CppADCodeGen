#ifndef CPPAD_CG_CUSTOM_POSITION_INCLUDED
#define CPPAD_CG_CUSTOM_POSITION_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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

/**
 * Useful class for storing sparsity patterns
 */
class CustomPosition {
private:
    bool filterDefined_;
    /// allowed elements
    std::vector<std::vector<bool> > elFilter_;
    bool fullDefined_;
    sparse_rc<CppAD::vector<size_t>> elements_;
public:

    inline CustomPosition() :
        filterDefined_(false),
        fullDefined_(false) {
    }

    template<class VectorSize>
    inline CustomPosition(size_t m, size_t n,
                          const VectorSize& filterRows,
                          const VectorSize& filterCols) :
        filterDefined_(true),
        elFilter_(m, std::vector<bool>(n, false)),
        fullDefined_(false) {
        CPPADCG_ASSERT_KNOWN(filterRows.size() == filterCols.size(), "The number of row indexes must be the same as the number of column indexes.")
        for (size_t i = 0; i < filterRows.size(); i++) {
            elFilter_[filterRows[i]][filterCols[i]] = true;
        }
    }

    template<class VectorSet>
    inline CustomPosition(size_t m, size_t n,
                          const VectorSet& filter) :
        filterDefined_(true),
        elFilter_(m, std::vector<bool>(n, false)),
        fullDefined_(false) {
        CPPADCG_ASSERT_KNOWN(filter.size() <= m, "Invalid number of rows.")

        for (size_t i = 0; i < filter.size(); i++) {
            for (size_t it : filter[i]) {
                elFilter_[i][it] = true;
            }
        }
    }

    inline bool isFilterDefined() const {
        return filterDefined_;
    }

    inline bool isFullDefined() const {
        return fullDefined_;
    }

    inline void setFullElements(const sparse_rc<CppAD::vector<size_t>>& elements) {
        elements_ = elements;
        filter(elements_);
        fullDefined_ = true;
    }

    inline const sparse_rc<CppAD::vector<size_t>>& getFullElements()const {
        return elements_;
    }

    inline void filter(sparse_rc<CppAD::vector<size_t>>& sparsity) const {
        if (!filterDefined_)
            return; // nothing to do

        size_t nnz = sparsity.nnz();
        auto& row = sparsity.row();
        auto& col = sparsity.col();

        sparse_rc<CppAD::vector<size_t>> filtered(sparsity.nr(), sparsity.nc(), nnz);

        size_t e = 0;
        for (size_t k = 0; k < nnz; ++k) {
            size_t i = row[k];
            size_t j = col[k];
            if (elFilter_[i][j]) {
                filtered.set(e, i,j);
                e++;
            }
        }

        if (nnz != e) {
            filtered.resize(sparsity.nr(), sparsity.nc(), e);
            sparsity = filtered;
        }
    }

};

} // END cg namespace
} // END CppAD namespace

#endif