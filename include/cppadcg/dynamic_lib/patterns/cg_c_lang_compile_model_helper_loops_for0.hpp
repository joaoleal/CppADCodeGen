#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_FOR0_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_FOR0_INCLUDED
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

    /***************************************************************************
     *  Methods related with loop insertion into the operation graph
     **************************************************************************/

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareForward0WithLoops(CodeHandler<Base>& handler,
                                                                 std::vector<CGBase>& y) {
        using namespace std;
        using CppAD::vector;

        std::vector<IndexedDependentLoopInfo<Base> > garbageCollection(y.size()); ////////// <- rethink this!!!!!!!
        garbageCollection.resize(y.size());

        map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > > dependentIndexes;
        // loop -> loop atomic evaluation -> results
        map<LoopAtomicFun<Base>*, map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> > > evaluations1it;

        for (size_t i = 0; i < y.size(); i++) {
            OperationNode<Base>* node = y[i].getOperationNode();
            if (node == NULL || node->getOperationType() != CGLoopAtomicResultOp) {
                continue;
            }

            OperationNode<Base>* loopEvalNode = node->getArguments()[0].getOperation();
            assert(loopEvalNode->getOperationType() == CGLoopForwardOp);

            LoopForwardOperationNode<Base>* loopForNode = static_cast<LoopForwardOperationNode<Base>*> (loopEvalNode);
            assert(loopForNode->getOrder() == 0);

            LoopAtomicFun<Base>* loop = &loopForNode->getLoopAtomicFun();

            const LoopPosition& depPos = loop->getTapeDependentIndex(i);
            size_t mTape = loop->getTapeDependentCount();
            size_t iteration = depPos.atomic / mTape;

            if (iteration == 0) {
                vector<OperationNode<Base>*>& allLoopEvalResults = evaluations1it[loop][loopForNode];
                allLoopEvalResults.resize(loop->getTapeDependentCount());
                allLoopEvalResults[depPos.tape] = node;
            }

            vector<IndexedDependentLoopInfo<Base>* >& depInfo = dependentIndexes[loop];
            depInfo.resize(mTape);
            IndexedDependentLoopInfo<Base>* origEl = depInfo[depPos.tape];
            if (origEl == NULL) {
                garbageCollection.resize(garbageCollection.size() + 1);
                origEl = &garbageCollection.back();
                depInfo[depPos.tape] = origEl;
            }
            origEl->indexes.resize(loop->getIterationCount());
            origEl->origVals.resize(loop->getIterationCount());
            origEl->indexes[iteration] = i;
            origEl->origVals[iteration] = y[i];
            origEl->pattern = loop->getDependentIndexPatterns()[depPos.tape];
        }

        prepareLoops(handler, y, evaluations1it, dependentIndexes);
    }

}

#endif
