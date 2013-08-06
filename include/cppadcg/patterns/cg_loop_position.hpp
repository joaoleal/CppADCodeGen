#ifndef CPPAD_CG_LOOP_POSITION_INCLUDED
#define CPPAD_CG_LOOP_POSITION_INCLUDED
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
     * Independent variable positions
     */
    class LoopPosition {
    public:
        size_t tape;
        size_t atomic;
        size_t original;

        inline LoopPosition() :
            tape(-1),
            atomic(-1),
            original(-1) {
        }

        inline LoopPosition(size_t t, size_t a, size_t o) :
            tape(t),
            atomic(a),
            original(o) {
        }
    };

    /**
     * Temporary variable positions
     */
    class LoopPositionTmp {
    public:
        size_t tape;
        size_t atomic;
        // the independent variables that this temporary variable depends on
        std::set<size_t> originalIndeps;

        inline LoopPositionTmp() :
            tape(-1),
            atomic(-1) {
        }

        inline LoopPositionTmp(size_t t, size_t a, const std::set<size_t>& o) :
            tape(t),
            atomic(a),
            originalIndeps(o) {
        }
    };

}

#endif