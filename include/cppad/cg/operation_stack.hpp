#ifndef CPPAD_CG_OPERATION_STACK_HPP
#define CPPAD_CG_OPERATION_STACK_HPP
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2019 Joao Leal
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
 * The steps for a node during the navigation through the graph of operation nodes
 */
enum class StackNavigationStep {
    Analyze, // analyze the current node (before visiting children)
    ChildrenVisited, // post process node after children have been visited
    Exit // leave node (go to the parent node)
};

/**
 * A stack element in the Operation Node stack.
 *
 * @tparam Base
 */
template<class Base>
class OperationStackData {
public:
    OperationNode<Base>* const node;
    size_t parentNodeScope;
    StackNavigationStep nextStep;

    inline OperationStackData(OperationNode<Base>& node,
                              size_t parentNodeScope) :
            node(&node),
            parentNodeScope(parentNodeScope),
            nextStep(StackNavigationStep::Analyze) {
    }
};

/**
 * A Stack of Operation Nodes.
 *
 * @tparam Base
 */
template<class Base>
class OperationStack {
private:
    std::vector<OperationStackData<Base>> _stack;
public:
    inline OperationStack() {
        _stack.reserve(100);
    }

    inline bool empty() const {
        return _stack.empty();
    }

    inline size_t size() const {
        return _stack.size();
    }

    inline void pop_back() {
        _stack.pop_back();
    }

    inline OperationStackData<Base>& back() {
        return _stack.back();
    }

    inline void emplace_back(OperationNode<Base>& node,
                             size_t parentNodeScope) {
        if (_stack.size() == _stack.capacity()) {
            _stack.reserve((_stack.size() * 3) / 2 + 1);
        }
        _stack.emplace_back(node, parentNodeScope);
    }

    inline OperationStackData<Base>& operator[](size_t i) {
        return _stack[i];
    }
};


} // END cg namespace
} // END CppAD namespace

#endif
