#ifndef CPPAD_CG_INDEX_OPERATION_NODE_INCLUDED
#define CPPAD_CG_INDEX_OPERATION_NODE_INCLUDED
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
     * An index reference operation node
     * 
     * @author Joao Leal
     */
    template<class Base>
    class IndexOperationNode : public OperationNode<Base> {
    private:
        const Index& index_;
    public:

        inline IndexOperationNode(const Index& index) :
            OperationNode<Base>(CGIndexOp),
            index_(index) {
        }

        inline IndexOperationNode(const Index& index,
                                  OperationNode<Base>& loopStartOrIndexAssign) :
            OperationNode<Base>(CGIndexOp, Argument<Base>(loopStartOrIndexAssign)),
            index_(index) {
            CPPADCG_ASSERT_KNOWN(loopStartOrIndexAssign.getOperationType() == CGLoopStartOp ||
                                 (loopStartOrIndexAssign.getOperationType() == CGIndexAssignOp &&
                                 &static_cast<IndexAssignOperationNode<Base>&> (loopStartOrIndexAssign).getIndex() == &index),
                                 "Invalid argument");
        }

        inline const Index& getIndex() const {
            return index_;
        }

        inline virtual ~IndexOperationNode() {
        }

    };

}

#endif