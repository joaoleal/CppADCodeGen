#ifndef CPPAD_CG_LOOP_FORWARD_OPERATION_NODE_INCLUDED
#define CPPAD_CG_LOOP_FORWARD_OPERATION_NODE_INCLUDED
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
     * A forward mode operation node for the evaluation of the a model in a loop.
     * 
     * This is a custom OperationNode class and therefore cannot be transformed
     * into any other node type (makeAlias() and setOperation() might not work).
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LoopForwardOperationNode : public LoopEvaluationOperationNode<Base> {
    public:

        inline LoopForwardOperationNode(LoopAtomicFun<Base>& loopFunc,
                                        size_t p, // order
                                        size_t resultSize, // dimension of results of this atomic function call
                                        size_t tapeResultSize, // dimension of results relative to the tape
                                        const std::vector<Argument<Base> >& tx) :
            LoopEvaluationOperationNode<Base>(CGLoopForwardOp, loopFunc, p, resultSize, tapeResultSize, tx) {
        }

        inline virtual ~LoopForwardOperationNode() {
        }

    };

}

#endif