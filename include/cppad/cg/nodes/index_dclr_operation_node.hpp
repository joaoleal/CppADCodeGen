#ifndef CPPAD_CG_INDEX_DCLR_OPERATION_NODE_INCLUDED
#define CPPAD_CG_INDEX_DCLR_OPERATION_NODE_INCLUDED
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
 * An index declaration operation node
 * 
 * This is a custom OperationNode class and therefore cannot be transformed
 * into any other node type (makeAlias() and setOperation() might not work).
 * 
 * @author Joao Leal
 */
template<class Base>
class IndexDclrOperationNode : public OperationNode<Base> {
public:

    inline IndexDclrOperationNode(const std::string& name) :
        OperationNode<Base>(CGOpCode::IndexDeclaration) {
        CPPADCG_ASSERT_KNOWN(!name.empty(), "index name cannot be empty");
        this->setName(name);
    }

    inline virtual ~IndexDclrOperationNode() {
    }

};

} // END cg namespace
} // END CppAD namespace

#endif