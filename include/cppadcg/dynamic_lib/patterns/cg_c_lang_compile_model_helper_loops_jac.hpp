#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_JAC_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_JAC_INCLUDED
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
    void CLangCompileModelHelper<Base>::prepareSparseJacobianWithLoops(CodeHandler<Base>& handler,
                                                                       std::vector<CGBase>& jac) {
        using namespace std;
        using CppAD::vector;

        size_t nnz = _jacSparsity.rows.size();

        std::vector<IndexedDependentLoopInfo<Base> > garbageCollection; ////////// <- rethink this!!!!!!!
        garbageCollection.reserve(nnz);

        /**
         * Generate index patterns for the jacobian elements resulting from loops
         */
        // loop -> equation pattern -> tape independent -> orig independent(temporaries only) -> iteration = position
        map<LoopAtomicFun<Base>*, map<size_t, map<TapeVarType, IndexedDependentLoopInfo<Base>* > > > jacIndexPatterns;

        map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > > dependentIndexes;
        // loop -> loop atomic evaluation -> results
        map<LoopAtomicFun<Base>*, map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> > > evaluations1it;


        for (size_t e = 0; e < nnz; e++) {
            size_t i = _jacSparsity.rows[e];
            size_t j = _jacSparsity.cols[e];
            CGBase& jacVal = jac[e];

            // find LOOP + get loop results
            LoopAtomicFun<Base>* loop = NULL;

            // loop evaluation -> results
            map<LoopEvaluationOperationNode<Base>*, map<size_t, OperationNode<Base>*> > evals;
            findLoopEvaluations(handler, jacVal.getOperationNode(), evals);

            if (evals.empty()) {
                // no choice but to check all loops
                const std::map<size_t, LoopAtomicFun<Base>*>& loops = handler.getLoops();
                typename std::map<size_t, LoopAtomicFun<Base>*>::const_iterator itl;
                for (itl = loops.begin(); itl != loops.end(); ++itl) {
                    LoopAtomicFun<Base>* l = itl->second;
                    const std::map<size_t, LoopPosition*>& depIndexes = l->getOriginalDependentIndexes();
                    std::map<size_t, LoopPosition*>::const_iterator iti = depIndexes.find(i);
                    if (iti != depIndexes.end()) {
                        loop = l;
                        break;
                    }
                }

                if (loop == NULL)
                    continue; // result in jacobian is not related with a loop
            } else {
                assert(evals.size() == 1);
                LoopEvaluationOperationNode<Base>* loopEvalNode = evals.begin()->first;

                if (loopEvalNode->getOperationType() == CGLoopForwardOp) {
                    assert(loopEvalNode->getOrder() == 1);
                } else {
                    assert(loopEvalNode->getOrder() == 0);
                }

                loop = &loopEvalNode->getLoopAtomicFun();
            }

            map<size_t, map<TapeVarType, IndexedDependentLoopInfo<Base>* > >& lRawIndexes = jacIndexPatterns[loop];

            const LoopPosition& depPos = loop->getTapeDependentIndex(i);
            size_t tapeI = depPos.tape;
            if (!loop->isIndependentVariablesFullyDetached(tapeI)) {
                throw CGException("There are independent variables which appear as indexed and non-indexed for the same equation pattern");
            }
            size_t iteration = depPos.atomic / loop->getTapeDependentCount();

            set<size_t> tapeJs = loop->getIndependentTapeIndexes(tapeI, j);
            assert(tapeJs.size() >= 1); // if >1 then the independent is used by several temporary variables

            bool isTemporary = loop->isTemporary(*tapeJs.begin());
            size_t jRef = isTemporary ? j : 0;

            TapeVarType tapeJ(*tapeJs.begin(), jRef);

            IndexedDependentLoopInfo<Base>* origJacEl = NULL;
            map<TapeVarType, IndexedDependentLoopInfo<Base>* >& ref = lRawIndexes[tapeI];
            typename map<TapeVarType, IndexedDependentLoopInfo<Base>* >::const_iterator it;
            it = ref.find(tapeJ);
            if (it != ref.end()) {
                origJacEl = it->second;
            }

            if (origJacEl == NULL) {
                // the vector will never allocate more space so this is safe:
                garbageCollection.resize(garbageCollection.size() + 1);
                origJacEl = &garbageCollection.back();
                origJacEl->indexes.resize(loop->getIterationCount(), nnz);
                origJacEl->origVals.resize(loop->getIterationCount());
                ref[tapeJ] = origJacEl;
                dependentIndexes[loop].push_back(origJacEl);
            }
            origJacEl->indexes[iteration] = e;
            origJacEl->origVals[iteration] = jacVal;

            if (iteration == 0) {
                if (evals.size() > 0) {
                    typename map<LoopEvaluationOperationNode<Base>*, map<size_t, OperationNode<Base>*> >::const_iterator itE;
                    for (itE = evals.begin(); itE != evals.end(); ++itE) {
                        LoopEvaluationOperationNode<Base>* loopEvalNode = itE->first;
                        const map<size_t, OperationNode<Base>*>& results = itE->second;

                        vector<OperationNode<Base>*>& allLoopEvalResults = evaluations1it[loop][loopEvalNode];

                        size_t result_size = loopEvalNode->getTapeResultSize(); // tape result size
                        allLoopEvalResults.resize(result_size);

                        typename map<size_t, OperationNode<Base>*>::const_iterator itR;
                        for (itR = results.begin(); itR != results.end(); ++itR) {
                            size_t tapeIndex = itR->second->getInfo()[1];
                            allLoopEvalResults[tapeIndex] = itR->second;
                        }
                    }
                }
            }
        }

        /**
         * Generate index patterns for the dependent variables
         */
        // loop loops :)
        typename map<LoopAtomicFun<Base>*, map<size_t, map<TapeVarType, IndexedDependentLoopInfo<Base>* > > >::iterator itl;
        for (itl = jacIndexPatterns.begin(); itl != jacIndexPatterns.end(); ++itl) {

            // loop equation patterns
            typename map<size_t, map<TapeVarType, IndexedDependentLoopInfo<Base>* > >::iterator itI;
            for (itI = itl->second.begin(); itI != itl->second.end(); ++itI) {

                // loop tape variables
                typename map<TapeVarType, IndexedDependentLoopInfo<Base>* >::iterator itJ;
                for (itJ = itI->second.begin(); itJ != itI->second.end(); ++itJ) {
                    IndexedDependentLoopInfo<Base>* jacEleInfo = itJ->second;

                    // make sure all elements are requested for all iterations
                    std::vector<size_t>::const_iterator ite;
                    for (ite = jacEleInfo->indexes.begin(); ite != jacEleInfo->indexes.end(); ++ite) {
                        if (*ite == nnz) {
                            throw CGException("All jacobian elements of an equation pattern (equation in a loop) must be requested for all iterations");
                        }
                    }

                    jacEleInfo->pattern = IndexPattern::detect(LoopAtomicFun<Base>::ITERATION_INDEX, jacEleInfo->indexes);
                    handler.manageLoopDependentIndexPattern(jacEleInfo->pattern);

                }
            }
        }

        prepareLoops(handler, jac, evaluations1it, dependentIndexes);
    }

}

#endif