#ifndef CPPAD_CG_LINEAR_INDEX_PATTERN_INCLUDED
#define CPPAD_CG_LINEAR_INDEX_PATTERN_INCLUDED
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
     * Linear pattern y = a * x + b
     */
    class LinearIndexPattern : public IndexPattern {
    protected:
        // slope
        long a_;
        // constant term
        long b_;
    public:

        inline LinearIndexPattern(const Index& index, long a, long b) :
            IndexPattern(index),
            a_(a),
            b_(b) {
        }

        inline long getLinearSlope() const {
            return a_;
        }

        inline long getLinearConstantTerm() const {
            return b_;
        }

        inline virtual IndexPatternType getType() const {
            return LINEAR;
        }

        inline virtual ~LinearIndexPattern() {
        }
    };

}

#endif