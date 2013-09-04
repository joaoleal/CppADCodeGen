#ifndef CPPAD_CG_CUSTOM_LOOP_NODE_INFO_INCLUDED
#define CPPAD_CG_CUSTOM_LOOP_NODE_INFO_INCLUDED
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
    class CustomLoopNodeInfo : public LoopNodeInfo<Base> {
    protected:
        const Index& index_;
        /**
         * Number of loop iterations
         */
        const size_t iterationCount_;
        IndexOperationNode<Base>* iterationCountNode_; // CGIndexOp

    public:

        inline CustomLoopNodeInfo(const Index& index,
                                  IndexOperationNode<Base>& iterationCountNode) :
            index_(index),
            iterationCount_(0),
            iterationCountNode_(&iterationCountNode) {
        }

        inline CustomLoopNodeInfo(const Index& index,
                                  size_t iterationCount) :
            index_(index),
            iterationCount_(iterationCount),
            iterationCountNode_(NULL) {
        }

        virtual const std::string* getLoopName() const {
            return NULL;
        }

        virtual inline const Index& getIndex() const {
            return index_;
        }

        virtual inline IndexOperationNode<Base>* getIterationCountNode() const {
            return iterationCountNode_;
        }

        virtual inline const size_t getIterationCount() const {
            return iterationCount_;
        }

        inline virtual ~CustomLoopNodeInfo() {
        }
    };
}

#endif