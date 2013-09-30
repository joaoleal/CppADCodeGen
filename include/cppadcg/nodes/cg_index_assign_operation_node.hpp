#ifndef CPPAD_CG_INDEX_ASSIGN_OPERATION_NODE_INCLUDED
#define CPPAD_CG_INDEX_ASSIGN_OPERATION_NODE_INCLUDED
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
     * An index reference operation node.
     * 
     * This is a custom OperationNode class and therefore cannot be transformed
     * into any other node type (makeAlias() and setOperation() might not work).
     * 
     * @author Joao Leal
     */
    template<class Base>
    class IndexAssignOperationNode : public OperationNode<Base> {
    private:
        const IndexPattern& indexPattern_;
    public:

        /**
         * 
         * @param index the index that is being assigned
         * @param indexPattern the index expression used to define the index value
         * @param index1
         */
        inline IndexAssignOperationNode(IndexDclrOperationNode<Base>& index,
                                        const IndexPattern& indexPattern,
                                        IndexOperationNode<Base>& index1) :
            OperationNode<Base>(CGIndexAssignOp, Argument<Base>(index), Argument<Base>(index1)),
            indexPattern_(indexPattern) {
        }

        inline IndexAssignOperationNode(IndexDclrOperationNode<Base>& index,
                                        const IndexPattern& indexPattern,
                                        IndexOperationNode<Base>* index1,
                                        IndexOperationNode<Base>* index2) :
            OperationNode<Base>(CGIndexAssignOp, std::vector<size_t> (0), createArguments(index, index1, index2)),
            indexPattern_(indexPattern) {
        }

        inline IndexDclrOperationNode<Base>& getIndex() const {
            const std::vector<Argument<Base> >& args = this->getArguments();
            CPPADCG_ASSERT_KNOWN(!args.empty(), "Invalid number of arguments");

            OperationNode<Base>* aNode = args[0].getOperation();
            CPPADCG_ASSERT_KNOWN(aNode != NULL && aNode->getOperationType() == CGIndexDeclarationOp, "Invalid argument operation type");

            return static_cast<IndexDclrOperationNode<Base>&> (*aNode);
        }

        inline const IndexPattern& getIndexPattern() const {
            return indexPattern_;
        }

        inline std::vector<const IndexDclrOperationNode<Base>*> getIndexPatternIndexes() const {
            std::vector<const IndexDclrOperationNode<Base>*> iargs;

            const std::vector<Argument<Base> >& args = this->getArguments();

            CPPADCG_ASSERT_KNOWN(args[1].getOperation() != NULL &&
                                 args[1].getOperation()->getOperationType() == CGIndexOp, "Invalid argument operation type");
            iargs.push_back(&static_cast<IndexOperationNode<Base>*> (args[1].getOperation())->getIndex());

            if (args.size() > 2) {
                CPPADCG_ASSERT_KNOWN(args[2].getOperation() != NULL &&
                                     args[2].getOperation()->getOperationType() == CGIndexOp, "Invalid argument operation type");
                iargs.push_back(&static_cast<IndexOperationNode<Base>*> (args[2].getOperation())->getIndex());
            }

            return iargs;
        }

        inline virtual ~IndexAssignOperationNode() {
        }

    private:

        inline static std::vector<Argument<Base> > createArguments(IndexDclrOperationNode<Base>& index,
                                                                   IndexOperationNode<Base>* index1,
                                                                   IndexOperationNode<Base>* index2) {
            std::vector<Argument<Base> > args(1 + (index1 != NULL)+ (index2 != NULL));

            args[0] = Argument<Base>(index);
            if (index1 != NULL)
                args[1] = Argument<Base>(*index1);
            if (index2 != NULL) {
                args.back() = Argument<Base>(*index2);
            }

            return args;
        }

    };

}
#endif