#ifndef CPPAD_CG_RANDOM_INDEX_PATTERN_INCLUDED
#define CPPAD_CG_RANDOM_INDEX_PATTERN_INCLUDED
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
     * Random pattern
     */
    class RandomIndexPattern : public IndexPattern {
    public:

        inline RandomIndexPattern() {
        }

        inline virtual IndexPatternType getType() const {
            return RANDOM;
        }

        inline virtual ~RandomIndexPattern() {
        }
    };

}

#endif