#ifndef CPPAD_CG_INDEX_PATTERN_IMPL_INCLUDED
#define CPPAD_CG_INDEX_PATTERN_IMPL_INCLUDED
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

namespace CppAD {

    template<class VectorSizeT>
    IndexPattern* IndexPattern::detect(const Index& indexX, const VectorSizeT& x2y) {
        assert(x2y.size() > 0);

        size_t maxCount = std::min(std::max(3ul, x2y.size() / 4), 8ul);
        std::map<size_t, IndexPattern*> linearSections = SectionedIndexPattern::detectLinearSections(indexX, x2y, maxCount);

        if (linearSections.size() == 1) {
            return linearSections.begin()->second;
        } else if (!linearSections.empty()) {
            return new SectionedIndexPattern(linearSections);
        } else {
            throw CGException("Random index patterns not implemented yet!");
            return new RandomIndexPattern(indexX);
        }

    }

    IndexPattern* IndexPattern::detect(const Index& indexX, const std::map<size_t, size_t>& x2y) {
        assert(!x2y.empty());

        size_t maxCount = std::min(std::max(3ul, x2y.size() / 4), 8ul);
        std::map<size_t, IndexPattern*> linearSections = SectionedIndexPattern::detectLinearSections(indexX, x2y, maxCount);

        if (linearSections.size() == 1) {
            return linearSections.begin()->second;
        } else if (!linearSections.empty()) {
            return new SectionedIndexPattern(linearSections);
        } else {
            throw CGException("Random index patterns not implemented yet!");
            return new RandomIndexPattern(indexX);
        }
    }

    inline bool IndexPattern::isConstant(const IndexPattern& ip) {
        if (ip.getType() == LINEAR) {
            const LinearIndexPattern& lip = static_cast<const LinearIndexPattern&> (ip);
            return lip.getLinearSlopeDy() == 0;
        }
        return false;
    }
}

#endif