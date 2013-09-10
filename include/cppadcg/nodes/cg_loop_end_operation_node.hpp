#ifndef CPPAD_CG_LOOP_END_OPERATION_NODE_INCLUDED
#define CPPAD_CG_LOOP_END_OPERATION_NODE_INCLUDED
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
     * An operation node that marks the end of a loop.
     * 
     * This is a custom OperationNode class and therefore cannot be transformed
     * into any other node type (makeAlias() and setOperation() might not work).
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LoopEndOperationNode : public OperationNode<Base> {
    private:
        LoopStartOperationNode<Base>& loopStart_;
        const LoopNodeInfo<Base>& loopInfo_;
    public:

        inline LoopEndOperationNode(const LoopNodeInfo<Base>& info,
                                    LoopStartOperationNode<Base>& loopStart,
                                    const std::vector<Argument<Base> >& endArgs) :
            OperationNode<Base>(CGLoopEndOp, std::vector<size_t>(0), createArguments(loopStart, endArgs)),
            loopStart_(loopStart),
            loopInfo_(info) {
        }

        inline const LoopStartOperationNode<Base>& getLoopStart() const {
#ifndef NDEBUG
            const std::vector<Argument<Base> >& args = this->getArguments();
            CPPADCG_ASSERT_KNOWN(args.size() > 0, "There must be at least one argument");
            CPPADCG_ASSERT_KNOWN(args[0].getOperation() == &loopStart_, "The first argument must be the loop start operation");
#endif

            return loopStart_;
        }

        inline const LoopNodeInfo<Base>& getLoopInfo() const {
            return loopInfo_;
        }

        inline virtual ~LoopEndOperationNode() {
        }

    private:

        static inline std::vector<Argument<Base> > createArguments(LoopStartOperationNode<Base>& lstart,
                                                                   const std::vector<Argument<Base> >& endArgs) {
            std::vector<Argument<Base> > args(1 + endArgs.size());
            args[0] = Argument<Base>(lstart);
            std::copy(endArgs.begin(), endArgs.end(), args.begin() + 1);
            return args;
        }

    };

}

#endif