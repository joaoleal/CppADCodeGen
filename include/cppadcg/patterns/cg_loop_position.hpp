#ifndef CPPAD_CG_LOOP_POSITION_INCLUDED
#define CPPAD_CG_LOOP_POSITION_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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

    /**
     * Variable positions
     */
    class LoopPosition {
    public:
        size_t tape;
        size_t original;

        inline LoopPosition() :
            tape(-1),
            original(-1) {
        }

        inline LoopPosition(size_t t, size_t o) :
            tape(t),
            original(o) {
        }
    };

    /**
     * Variable position which changes from iteration to iteration
     */
    class LoopIndexedPosition : public LoopPosition {
    public:
        size_t iteration;

        inline LoopIndexedPosition() :
            LoopPosition(),
            iteration(-1) {
        }

        inline LoopIndexedPosition(size_t t, size_t o, size_t it) :
            LoopPosition(t, o),
            iteration(it) {
        }
    };
}

#endif