#ifndef CPPAD_CG_LOOP_EVALUATION_OPERATION_NODE_INCLUDED
#define CPPAD_CG_LOOP_EVALUATION_OPERATION_NODE_INCLUDED
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
     * An operation node for the evaluation of the a model in a loop.
     * 
     * This is a custom OperationNode class and therefore cannot be transformed
     * into any other node type (makeAlias() and setOperation() might not work).
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LoopEvaluationOperationNode : public OperationNode<Base> {
    private:
        LoopAtomicFun<Base>& loopFunc_;
        const size_t p_;
        const size_t resultSize_; // dimension of results of this atomic function call
        const size_t tapeResultSize_; // dimension of results relative to the tape
    public:

        inline LoopAtomicFun<Base>& getLoopAtomicFun() {
            return loopFunc_;
        }

        inline const LoopAtomicFun<Base>& getLoopAtomicFun() const {
            return loopFunc_;
        }

        inline size_t getOrder() const {
            return p_;
        }

        inline size_t getResultSize() const {
            return resultSize_;
        }

        inline size_t getTapeResultSize() const {
            return tapeResultSize_;
        }

        inline virtual ~LoopEvaluationOperationNode() {
        }

    protected:

        inline LoopEvaluationOperationNode(CGOpCode op,
                                           LoopAtomicFun<Base>& loopFunc,
                                           size_t p, // order
                                           size_t resultSize, // dimension of results of this atomic function call
                                           size_t tapeResultSize, // dimension of results relative to the tape
                                           const std::vector<Argument<Base> >& args) :
            OperationNode<Base>(op, std::vector<size_t>(0), args),
            loopFunc_(loopFunc),
            p_(p),
            resultSize_(resultSize),
            tapeResultSize_(tapeResultSize) {
            assert(op == CGLoopForwardOp || op == CGLoopReverseOp);
        }

    };

}

#endif