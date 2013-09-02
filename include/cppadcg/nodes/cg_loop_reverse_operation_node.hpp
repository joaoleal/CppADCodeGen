#ifndef CPPAD_CG_LOOP_REVERSE_OPERATION_NODE_INCLUDED
#define CPPAD_CG_LOOP_REVERSE_OPERATION_NODE_INCLUDED
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
     * A loop end operation node
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LoopReverseOperationNode : public LoopEvaluationOperationNode<Base> {
    public:

        inline LoopReverseOperationNode(LoopAtomicFun<Base>& loopFunc,
                                        size_t p, // order
                                        size_t resultSize, // dimension of results of this atomic function call
                                        size_t tapeResultSize, // dimension of results relative to the tape
                                        const std::vector<Argument<Base> >& tx,
                                        const std::vector<Argument<Base> >& py) :
            LoopEvaluationOperationNode<Base>(CGLoopReverseOp, loopFunc, p, resultSize, tapeResultSize, combine(tx, py)) {
        }

        inline virtual ~LoopReverseOperationNode() {
        }

    private:

        static inline std::vector<Argument<Base> > combine(const std::vector<Argument<Base> >& tx,
                                                           const std::vector<Argument<Base> >& py) {
            std::vector<Argument<Base> > args(tx.size() + py.size());
            std::copy(tx.begin(), tx.end(), args.begin());
            std::copy(py.begin(), py.end(), args.begin() + tx.size());
            return args;
        }
    };

}

#endif