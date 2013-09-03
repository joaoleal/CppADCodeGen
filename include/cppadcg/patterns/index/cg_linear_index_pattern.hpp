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
     * Linear pattern y = ((x - offset) / dx) * dy + b
     */
    class LinearIndexPattern : public IndexPattern {
    protected:
        long xOffset_;
        // slope
        long dy_;
        long dx_;
        // constant term
        long b_;
    public:

        inline LinearIndexPattern(const Index& index, long xOffset, long dy, long dx, long b) :
            IndexPattern(index),
            xOffset_(xOffset),
            dy_(dy),
            dx_(dx),
            b_(b) {
        }

        inline long getXOffset()const {
            return xOffset_;
        }

        inline long getLinearSlopeDy() const {
            return dy_;
        }

        inline long getLinearSlopeDx() const {
            return dx_;
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