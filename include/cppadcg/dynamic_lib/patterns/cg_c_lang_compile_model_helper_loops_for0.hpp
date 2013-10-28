#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_FOR0_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_FOR0_INCLUDED
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

    /***************************************************************************
     *  Methods related with loop insertion into the operation graph
     **************************************************************************/

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::prepareForward0WithLoops(CodeHandler<Base>& handler,
                                                                              const vector<CGBase>& x) {
        using namespace std;
        using namespace loops;
        using CppAD::vector;

        vector<CGBase> y(_fun.Range());

        // temporaries
        vector<CGBase> tmps;

        /**
         * original equations outside the loops 
         */
        if (_funNoLoops != NULL) {
            const std::vector<size_t>& origEq = _funNoLoops->getOrigDependentIndexes();

            vector<CGBase> depNL = _funNoLoops->getTape().Forward(0, x);

            // original equations
            for (size_t e = 0; e < origEq.size(); e++) {
                y[origEq[e]] = depNL[e];
            }

            tmps.resize(depNL.size() - origEq.size());
            for (size_t i = origEq.size(); i < depNL.size(); i++)
                tmps[i - origEq.size()] = depNL[i];
        }

        /**
         * equations in loops
         */
        IndexDclrOperationNode<Base>* iterationIndexDcl = new IndexDclrOperationNode<Base>(LoopModel<Base>::ITERATION_INDEX_NAME);
        handler.manageOperationNodeMemory(iterationIndexDcl);

        typename std::set<LoopModel<Base>* >::const_iterator itl;
        size_t l = 0;
        for (itl = _loopTapes.begin(); itl != _loopTapes.end(); ++itl, l++) {
            LoopModel<Base>& lModel = **itl;
            size_t nIterations = lModel.getIterationCount();
            const std::vector<std::vector<LoopPosition> >& dependents = lModel.getDependentIndexes();

            /**
             * make the loop start
             */
            LoopStartOperationNode<Base>* loopStart = new LoopStartOperationNode<Base>(*iterationIndexDcl, nIterations);
            handler.manageOperationNodeMemory(loopStart);

            IndexOperationNode<Base>* iterationIndexOp = new IndexOperationNode<Base>(*loopStart);
            handler.manageOperationNodeMemory(iterationIndexOp);
            std::set<IndexOperationNode<Base>*> indexesOps;
            indexesOps.insert(iterationIndexOp);

            /**
             * evaluate the loop body
             */
            vector<CGBase> indexedIndeps = createIndexedIndependents(handler, lModel, *iterationIndexOp);
            vector<CGBase> xl = createLoopIndependentVector(handler, lModel, indexedIndeps, x, tmps);
            if (xl.size() == 0) {
                xl.resize(1); // does not depend on any variable but CppAD requires at least one
                xl[0] = Base(0);
            }
            vector<CGBase> yl = lModel.getTape().Forward(0, xl);

            /**
             * make the loop end
             */
            const vector<IndexPattern*>& depPatterns = lModel.getDependentIndexPatterns();
            vector<std::pair<CGBase, IndexPattern*> > indexedLoopResults(yl.size());
            for (size_t i = 0; i < yl.size(); i++) {
                indexedLoopResults[i] = std::make_pair(yl[i], depPatterns[i]);
            }
            size_t assignOrAdd = 0;
            LoopEndOperationNode<Base>* loopEnd = createLoopEnd(handler, *loopStart, indexedLoopResults, indexesOps, assignOrAdd);

            std::vector<size_t> info(1);
            std::vector<Argument<Base> > args(1);
            for (size_t i = 0; i < dependents.size(); i++) {
                for (size_t it = 0; it < nIterations; it++) {
                    // an additional alias variable is required so that each dependent variable can have its own ID
                    size_t e = dependents[i][it].original;
                    info[0] = e;
                    args[0] = Argument<Base>(*loopEnd);
                    y[e] = handler.createCG(new OperationNode<Base> (CGDependentRefRhsOp, info, args));
                }
            }

            /**
             * move non-indexed expressions outside loop
             */
            moveNonIndexedOutsideLoop(*loopStart, *loopEnd);
        }

        return y;
    }

}

#endif
