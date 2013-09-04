#ifndef CPPAD_CG_LOOP_START_OPERATION_NODE_INCLUDED
#define CPPAD_CG_LOOP_START_OPERATION_NODE_INCLUDED
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
     * A loop start operation node
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LoopStartOperationNode : public OperationNode<Base> {
    private:
        const LoopNodeInfo<Base>& loopInfo_;
    public:

        inline LoopStartOperationNode(const LoopNodeInfo<Base>& info) :
            OperationNode<Base>(CGLoopStartOp),
            loopInfo_(info) {
            if (info.getIterationCountNode() != NULL) {
                this->getArguments().push_back(Argument<Base>(*loopInfo_.getIterationCountNode()));
            }
        }

        inline const LoopNodeInfo<Base>& getLoopInfo() const {
            return loopInfo_;
        }

        inline virtual ~LoopStartOperationNode() {
        }

    };

}

#endif