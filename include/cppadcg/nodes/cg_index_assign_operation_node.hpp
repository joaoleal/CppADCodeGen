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
     * An index reference operation node
     * 
     * @author Joao Leal
     */
    template<class Base>
    class IndexAssignOperationNode : public OperationNode<Base> {
    private:
        const Index& index_;
        const IndexPattern& indexPattern_;
    public:

        inline IndexAssignOperationNode(const Index& index,
                                        const IndexPattern& indexPattern,
                                        IndexOperationNode<Base>& index1) :
            OperationNode<Base>(CGIndexAssignOp, Argument<Base>(index1)),
            index_(index),
            indexPattern_(indexPattern) {
        }

        inline IndexAssignOperationNode(const Index& index,
                                        const IndexPattern& indexPattern,
                                        IndexOperationNode<Base>* index1,
                                        IndexOperationNode<Base>* index2) :
            OperationNode<Base>(CGIndexAssignOp, std::vector<size_t> (0), createArguments(index1, index2)),
            index_(index),
            indexPattern_(indexPattern) {
        }

        inline const Index& getIndex() const {
            return index_;
        }

        inline const IndexPattern& getIndexPattern() const {
            return indexPattern_;
        }

        inline virtual ~IndexAssignOperationNode() {
        }

    private:

        inline static std::vector<Argument<Base> > createArguments(IndexOperationNode<Base>* index1,
                                                                   IndexOperationNode<Base>* index2) {
            std::vector<Argument<Base> > args((index1 != NULL)+ (index2 != NULL));

            if (index1 != NULL)
                args[0] = Argument<Base>(*index1);
            if (index2 != NULL) {
                args.back() = Argument<Base>(*index2);
            }

            return args;
        }

    };

}
#endif