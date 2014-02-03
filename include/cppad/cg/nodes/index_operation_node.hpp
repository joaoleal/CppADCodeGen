#ifndef CPPAD_CG_INDEX_OPERATION_NODE_INCLUDED
#define CPPAD_CG_INDEX_OPERATION_NODE_INCLUDED
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
namespace cg {

/**
 * An index reference operation node
 * 
 * This is a custom OperationNode class and therefore cannot be transformed
 * into any other node type (makeAlias() and setOperation() might not work).
 * 
 * @author Joao Leal
 */
template<class Base>
class IndexOperationNode : public OperationNode<Base> {
public:

    inline IndexOperationNode(IndexDclrOperationNode<Base>& indexDcl) :
        OperationNode<Base>(CGIndexOp, indexDcl) {
    }

    inline IndexOperationNode(LoopStartOperationNode<Base>& loopStart) :
        OperationNode<Base>(CGIndexOp,{loopStart.getIndex(), loopStart}) {
    }

    inline IndexOperationNode(IndexAssignOperationNode<Base>& indexAssign) :
        OperationNode<Base>(CGIndexOp,{indexAssign.getIndex(), indexAssign}) {
    }

    inline bool isDefinedLocally() const {
        return this->getArguments().size() > 1;
    }

    inline const IndexDclrOperationNode<Base>& getIndex() const {
        const std::vector<Argument<Base> >& args = this->getArguments();
        CPPADCG_ASSERT_KNOWN(!args.empty(), "Invalid number of arguments");

        OperationNode<Base>* aNode = args[0].getOperation();
        CPPADCG_ASSERT_KNOWN(aNode != nullptr && aNode->getOperationType() == CGIndexDeclarationOp, "Invalid argument operation type");

        return static_cast<const IndexDclrOperationNode<Base>&> (*aNode);
    }

    inline void makeAssigmentDependent(IndexAssignOperationNode<Base>& indexAssign) {
        std::vector<Argument<Base> >& args = this->getArguments();

        args.resize(2);
        args[0] = indexAssign.getIndex();
        args[1] = indexAssign;
    }

    inline virtual ~IndexOperationNode() {
    }

};

} // END cg namespace
} // END CppAD namespace

#endif