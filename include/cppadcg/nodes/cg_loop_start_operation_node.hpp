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
     * An operation node that marks the beginning of a loop.
     * 
     * This is a custom OperationNode class and therefore cannot be transformed
     * into any other node type (makeAlias() and setOperation() might not work).
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LoopStartOperationNode : public OperationNode<Base> {
    public:

        inline LoopStartOperationNode(IndexDclrOperationNode<Base>& indexDcl, size_t iterationCount) :
            OperationNode<Base>(CGLoopStartOp, Argument<Base>(indexDcl)) {
            this->getInfo().push_back(iterationCount);
        }

        inline LoopStartOperationNode(IndexDclrOperationNode<Base>& indexDcl, IndexOperationNode<Base>& iterCount) :
            OperationNode<Base>(CGLoopStartOp, Argument<Base>(indexDcl), Argument<Base>(iterCount)) {
        }

        inline IndexDclrOperationNode<Base>& getIndex() const {
            const std::vector<Argument<Base> >& args = this->getArguments();
            CPPADCG_ASSERT_KNOWN(!args.empty(), "Invalid number of arguments");

            OperationNode<Base>* aNode = args[0].getOperation();
            CPPADCG_ASSERT_KNOWN(aNode != NULL && aNode->getOperationType() == CGIndexDeclarationOp, "Invalid argument operation type");

            return static_cast<IndexDclrOperationNode<Base>&> (*aNode);
        }

        inline IndexOperationNode<Base>* getIterationCountNode() const {
            if (this->getInfo().empty()) {
                CPPADCG_ASSERT_KNOWN(this->getArguments().size() > 1, "Invalid number of arguments.");

                OperationNode<Base>* aNode = this->getArguments()[1].getOperation();
                CPPADCG_ASSERT_KNOWN(aNode != NULL && aNode->getOperationType() == CGIndexOp, "Invalid argument node type");

                return static_cast<IndexOperationNode<Base>*> (aNode);
            }

            return NULL;
        }

        inline const size_t getIterationCount() const {
            if (this->getInfo().empty()) {
                return 0;
            }
            return this->getInfo()[0];
        }

        inline virtual ~LoopStartOperationNode() {
        }

    };

}

#endif