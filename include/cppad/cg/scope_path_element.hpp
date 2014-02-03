#ifndef CPPAD_CG_SCOPE_PATH_INCLUDED
#define CPPAD_CG_SCOPE_PATH_INCLUDED
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

template<class Base>
class ScopePathElement {

public:
    // the color/index associated with the scope
    size_t color;
    // the node that marks the beginning of this scope
    OperationNode<Base>* beginning;
    // the node that marks the end of this scope
    OperationNode<Base>* end;
public:

    inline ScopePathElement(size_t color_ = 0, OperationNode<Base>* nEnd = nullptr, OperationNode<Base>* nBegin = nullptr) :
        color(color_),
        beginning(nBegin),
        end(nEnd) {
    }

};

} // END cg namespace
} // END CppAD namespace

#endif