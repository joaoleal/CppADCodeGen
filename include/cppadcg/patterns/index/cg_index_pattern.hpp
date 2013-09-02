#ifndef CPPAD_CG_INDEX_PATTERN_INCLUDED
#define CPPAD_CG_INDEX_PATTERN_INCLUDED
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

    /**
     * Generic index pattern
     */
    class IndexPattern {
    protected:
        std::vector<const Index*> indexes_;
    public:

        inline IndexPattern() {
        }

        inline IndexPattern(const Index& index) :
            indexes_(1) {
            indexes_[0] = &index;
        }

        virtual IndexPatternType getType() const = 0;

        inline const std::vector<const Index*>& getIndexes() const {
            return indexes_;
        }

        inline virtual ~IndexPattern() {
        }

        /***********************************************************************
         *   static methods
         **********************************************************************/
        /**
         * Detects the index pattern for the provided points
         * 
         * @param indexes maps the independents to the dependents (indexes[x] = y )
         * @return the generated index pattern (must be deleted by user)
         */
        template<class VectorSizeT>
        static inline IndexPattern* detect(const Index& index, const VectorSizeT& indexes);

        /**
         * Detects the index pattern for the provided points
         * 
         * @param indexes maps the independents to the dependents (x,y)
         * @return the generated index pattern (must be deleted by user)
         */
        static inline IndexPattern* detect(const Index& index, const std::map<size_t, size_t>& indexes);
    };

}

#endif