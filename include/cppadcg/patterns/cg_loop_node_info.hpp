#ifndef CPPAD_CG_LOOP_NODE_INFO_INCLUDED
#define CPPAD_CG_LOOP_NODE_INFO_INCLUDED
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
     * Information for a loop start operation node
     * 
     * @author Joao Leal
     */
    template <class Base>
    class LoopNodeInfo {
    private:
        const size_t loopId_;
    public:

        inline LoopNodeInfo() :
            loopId_(createNewLoopId()) {
        }

        /**
         * Provides a unique identifier for this loop.
         * 
         * @return a unique identifier ID
         */
        inline size_t getLoopId() const {
            return loopId_;
        }

        virtual const std::string* getLoopName() const = 0;

        virtual const Index& getIndex() const = 0;

        virtual IndexOperationNode<Base>* getIterationCountNode() const = 0;

        virtual const size_t getIterationCount() const = 0;

        inline virtual ~LoopNodeInfo() {
        }

    private:

        static size_t createNewLoopId() {
            CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL;
            static size_t count = 0;
            count++;
            return count;
        }
    };
}

#endif