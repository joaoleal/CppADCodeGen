#ifndef CPPAD_CG_SCOPE_PATH_INCLUDED
#define CPPAD_CG_SCOPE_PATH_INCLUDED
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

    template<class Base>
    class ScopePathElement {
    public:
        // the color/index associated with the scope
        size_t color;
        // the node that marks the beginning of this scope
        OperationNode<Base>* beginning;

        inline ScopePathElement(size_t color_ = 0, OperationNode<Base>* node = NULL) :
            color(color_),
            beginning(node) {
        }

    };

}

#endif
