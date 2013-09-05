#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_HESS_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_HESS_INCLUDED
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
    std::vector<CG<Base> > CLangCompileModelHelper<Base>::prepareSparseHessianWithLoops(CodeHandler<Base>& handler,
                                                                                        vector<CGBase>& indVars,
                                                                                        vector<CGBase>& w,
                                                                                        const std::vector<size_t>& lowerHessRows,
                                                                                        const std::vector<size_t>& lowerHessCols,
                                                                                        const std::vector<size_t>& lowerHessOrder,
                                                                                        const std::map<size_t, size_t>& duplicates) {
        size_t m = _fun->Range();
        //size_t n = _fun->Domain();

        std::vector<CGBase> hess(_hessSparsity.rows.size());

        size_t mLoopTotal = 0;
        std::vector<bool> eqLoop(m, false);
        typename std::set<LoopAtomicFun<Base>* >::const_iterator itl;
        for (itl = _loopAtomics.begin(); itl != _loopAtomics.end(); ++itl) {
            LoopAtomicFun<Base>* loop = *itl;

            size_t iterations = loop->getIterationCount();
            mLoopTotal += loop->getLoopDependentCount();

            const std::vector<std::vector<LoopPosition> >& depIndexes = loop->getDependentIndexes();
            for (size_t i = 0; i < depIndexes.size(); i++) {
                for (size_t iter = 0; iter < iterations; iter++) {
                    const LoopPosition& pos = depIndexes[i][iter];
                    eqLoop[pos.original] = true;
                }
            }
        }

        if (mLoopTotal < m) {
            /**
             * equations not in loops
             *  (must come before the loops because of the assigments to hess!)
             */
            vector<CGBase> ww(m);
            for (size_t i = 0; i < m; i++) {
                ww[i] = eqLoop[i] ? Base(0) : w[i];
            }

            CppAD::sparse_hessian_work work;
            vector<CGBase> lowerHess(lowerHessRows.size());
            _fun->SparseHessian(indVars, ww, _hessSparsity.sparsity, lowerHessRows, lowerHessCols, lowerHess, work);

            for (size_t i = 0; i < lowerHessOrder.size(); i++) {
                hess[lowerHessOrder[i]] = lowerHess[i];
            }

            // make use of the symmetry of the Hessian in order to reduce operations
            std::map<size_t, size_t>::const_iterator it2;
            for (it2 = duplicates.begin(); it2 != duplicates.end(); ++it2) {
                hess[it2->first] = hess[it2->second];
            }
        }

        /**
         * loops
         */
        std::map<LoopAtomicFun<Base>*, std::vector<CGBase> > loopHess; // hessian elements only for the equation in a given loop

        for (itl = _loopAtomics.begin(); itl != _loopAtomics.end(); ++itl) {
            LoopAtomicFun<Base>* loop = *itl;
            size_t iterations = loop->getIterationCount();

            vector<CGBase> wLoop(m);
            for (size_t i = 0; i < m; i++) {
                wLoop[i] = Base(0);
            }

            const std::vector<std::vector<LoopPosition> >& depIndexes = loop->getDependentIndexes();
            for (size_t i = 0; i < depIndexes.size(); i++) {
                for (size_t iter = 0; iter < iterations; iter++) {
                    const LoopPosition& pos = depIndexes[i][iter];
                    wLoop[pos.original] = w[pos.original];
                }
            }

            CppAD::sparse_hessian_work work;
            vector<CGBase> lowerHess(lowerHessRows.size());
            _fun->SparseHessian(indVars, wLoop, _hessSparsity.sparsity, lowerHessRows, lowerHessCols, lowerHess, work);

            std::vector<CGBase>& hessLoopl = loopHess[loop];
            hessLoopl.resize(_hessSparsity.rows.size());
            for (size_t e = 0; e < lowerHessOrder.size(); e++) {
                hessLoopl[lowerHessOrder[e]] = lowerHess[e];
            }

            // make use of the symmetry of the Hessian in order to reduce operations
            std::map<size_t, size_t>::const_iterator it2;
            for (it2 = duplicates.begin(); it2 != duplicates.end(); ++it2) {
                hessLoopl[it2->first] = hessLoopl[it2->second];
            }

            prepareSparseHessianForLoop(handler, loop, hessLoopl);

            for (size_t e = 0; e < hess.size(); e++) {
                hess[e] += hessLoopl[e];
            }
        }

        return hess;
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareSparseHessianForLoop(CodeHandler<Base>& handler,
                                                                    LoopAtomicFun<Base>* loop,
                                                                    std::vector<CGBase>& hess) {
        using namespace std;
        using CppAD::vector;

        size_t iterationCount = loop->getIterationCount();
        size_t nnz = _hessSparsity.rows.size();

        if (!loop->isTapeIndependentsFullyDetached2ndOrder()) {
            throw CGException("There are independent variables which appear as indexed and non-indexed for the same equation pattern");
        }

        std::vector<IndexedDependentLoopInfo<Base> > garbageCollection; ////////// <- rethink this!!!!!!!
        garbageCollection.reserve(nnz);

        /**
         * Generate index patterns for the hessian elements resulting from loops
         */
        // loop -> tape independent 1 -> orig independent(temporaries only) 1 ->
        //      -> tape independent 2 -> orig independent(temporaries only) 2 -> iteration = position
        map<TapeVarType, map<TapeVarType, IndexedDependentLoopInfo<Base>* > > hessIndexPatterns;

        map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > > dependentIndexes;
        // loop -> loop atomic evaluation -> results
        map<LoopAtomicFun<Base>*, map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> > > evaluations1it;

#if 0
        /**
         * Determine the Hessian pattern for the equations in the loop
         */
        vector<std::set<size_t> > loopHessPattern;
        {
            const std::vector<std::vector<LoopPosition> >& depIndexes = loop->getDependentIndexes();

            size_t n = _fun->Domain();

            /**
             * Determine the sparsity pattern p for Hessian of w^T F
             */
            vector<std::set<size_t> > r(n); // identity matrix
            for (size_t j = 0; j < n; j++)
                r[j].insert(j);
            _fun->ForSparseJac(n, r);

            vector<std::set<size_t> > s(1);
            for (size_t i = 0; i < depIndexes.size(); i++) {
                for (size_t it = 0; it < iterationCount; it++) {
                    s[0].insert(depIndexes[i][it].original);
                }
            }

            loopHessPattern = _fun->RevSparseHes(n, s, false);
        }

        printSparsityPattern(loopHessPattern, "hessian");
#endif

        for (size_t e = 0; e < nnz; e++) {
            CGBase& hessVal = hess[e];
            if (IdenticalZero(hessVal))
                continue;

            // loop evaluation -> result argument index -> results
            map<LoopEvaluationOperationNode<Base>*, map<size_t, OperationNode<Base>*> > evals;
            findLoopEvaluations(handler, loop, hessVal.getOperationNode(), evals); ////////////

            if (evals.empty())
                continue;

            if (evals.size() > 1 || evals.begin()->second.size() > 1) {
                throw CGException("Unable to generate expression for an hessian element which "
                                  "is either associated with multiple indexed independents "
                                  "or an indepedent with constant index.");
            }

#ifndef NDEBUG
            {
                LoopEvaluationOperationNode<Base>* loopEvalNode = evals.begin()->first;
                assert(loopEvalNode->getOperationType() == CGLoopReverseOp);
                assert(loopEvalNode->getOrder() == 1);
            }
#endif
            size_t j1 = _hessSparsity.rows[e];
            size_t j2 = _hessSparsity.cols[e];

            /**
             * If tapeJ1s.size()>1 then the independent is used by several
             * temporary variables. We know that it is NOT a mix of 
             * temporaries and indexed independents because this situation
             * is checked before.
             */
            set<size_t> tapeJ1s = loop->getIndependentTapeIndexes(j1);
            assert(tapeJ1s.size() >= 1);
            bool isTemporary1 = loop->isTemporary(*tapeJ1s.begin());
            size_t jRef1 = isTemporary1 ? j1 : 0;
            TapeVarType jTape1(*tapeJ1s.begin(), jRef1);

            set<size_t> tapeJ2s = loop->getIndependentTapeIndexes(j2);
            assert(tapeJ2s.size() >= 1); // if >1 then the independent is used by several temporary variables
            bool isTemporary2 = loop->isTemporary(*tapeJ2s.begin());
            size_t jRef2 = isTemporary2 ? j2 : 0;
            TapeVarType jTape2(*tapeJ2s.begin(), jRef2);

            IndexedDependentLoopInfo<Base>* origHessEl = NULL;

            map<TapeVarType, IndexedDependentLoopInfo<Base>* >& ref = hessIndexPatterns[jTape1];
            if (ref.size() != 0) {
                typename map<TapeVarType, IndexedDependentLoopInfo<Base>* >::const_iterator it;
                it = ref.find(jTape2);
                if (it != ref.end()) {
                    origHessEl = it->second;
                }
            }

            if (origHessEl == NULL) {
                // the vector will never allocate more space so this is safe:
                garbageCollection.resize(garbageCollection.size() + 1);
                origHessEl = &garbageCollection.back();
                origHessEl->indexes.resize(loop->getIterationCount(), nnz);
                origHessEl->origVals.resize(loop->getIterationCount());
                ref[jTape2] = origHessEl;
                dependentIndexes[loop].push_back(origHessEl);
            }

            /**
             * Determine the iteration
             */
            size_t iteration;

            bool indexed1 = loop->isIndexedIndependent(*tapeJ1s.begin());
            bool indexed2 = loop->isIndexedIndependent(*tapeJ2s.begin());
            if (indexed1 || indexed2) {
                std::set<size_t> iterations;
                if (indexed1) {
                    iterations = loop->getIterationsOfIndexedIndep(*tapeJ1s.begin(), j1);
                    assert(!iterations.empty());
                }

                if (!indexed1 || (iterations.size() > 1 && indexed2)) {
                    std::set<size_t> iterations2 = loop->getIterationsOfIndexedIndep(*tapeJ2s.begin(), j2);
                    if (iterations.empty()) {
                        iterations = iterations2;
                    } else {
                        std::set<size_t> intersection;
                        std::set_intersection(iterations.begin(), iterations.end(),
                                              iterations2.begin(), iterations2.end(),
                                              std::inserter(intersection, intersection.end()));

                        iterations.swap(intersection);
                    }
                }

                std::set<size_t>::const_iterator itIt;
                for (itIt = iterations.begin(); itIt != iterations.end(); ++itIt) {
                    size_t nIt = *itIt;
                    //size_t iteration = depPos.atomic / loop->getTapeDependentCount();

                    origHessEl->indexes[nIt] = e;
                    origHessEl->origVals[nIt] = hessVal;
                }

                if (iterations.size() == 1) {
                    iteration = *iterations.begin();
                } else {
                    if (iterations.find(0) != iterations.end()) {
                        iteration = 0;
                        // must change the argument of the call to the atomic
                        // function so that it only corresponds to a single iteration
                        throw CGException("Not implemented yet");
                    } else {
                        iteration = *iterations.begin();
                    }
                }

            } else {
                /**
                 * hessian element relative to two non-indexed independents
                 * (present in all iterations)
                 */
                iteration = 0; // it is actually present in all iterations!!!

                for (size_t iter = 0; iter < iterationCount; iter++) {
                    origHessEl->indexes[iter] = e;
                    origHessEl->origVals[iter] = hessVal;
                }
                // must change the argument of the call to the atomic
                // function so that it only corresponds to a single iteration
                throw CGException("Not implemented yet");
            }

            if (iteration == 0) {
                typename map<LoopEvaluationOperationNode<Base>*, map<size_t, OperationNode<Base>*> >::const_iterator itE;
                for (itE = evals.begin(); itE != evals.end(); ++itE) {
                    LoopEvaluationOperationNode<Base>* loopEvalNode = itE->first;
                    const map<size_t, OperationNode<Base>*>& results = itE->second;

                    assert(loopEvalNode->getOperationType() == CGLoopReverseOp);
                    size_t result_size = loopEvalNode->getTapeResultSize(); // tape result size

                    vector<OperationNode<Base>*>& allLoopEvalResults = evaluations1it[loop][loopEvalNode];
                    allLoopEvalResults.resize(result_size);

                    typename map<size_t, OperationNode<Base>*>::const_iterator itR;
                    for (itR = results.begin(); itR != results.end(); ++itR) {
                        size_t tapeIndex = itR->second->getInfo()[1]; ///////////////////////////////////
                        allLoopEvalResults[tapeIndex] = itR->second;
                    }
                }
            }
        }

        /**
         * Generate index patterns for the dependent variables
         */
        typename map<TapeVarType, map<TapeVarType, IndexedDependentLoopInfo<Base>* > >::iterator itJ1;
        for (itJ1 = hessIndexPatterns.begin(); itJ1 != hessIndexPatterns.end(); ++itJ1) {

            typename map<TapeVarType, IndexedDependentLoopInfo<Base>* >::iterator itJ2;
            for (itJ2 = itJ1->second.begin(); itJ2 != itJ1->second.end(); ++itJ2) {

                IndexedDependentLoopInfo<Base>* orig = itJ2->second;

                // make sure all element are requested
                std::vector<size_t>::const_iterator ite;
                for (ite = orig->indexes.begin(); ite != orig->indexes.end(); ++ite) {
                    if (*ite == nnz) {
                        throw CGException("All hessian elements must be requested for all iterations");
                    }
                }

                orig->pattern = IndexPattern::detect(LoopAtomicFun<Base>::ITERATION_INDEX, orig->indexes);
                handler.manageLoopDependentIndexPattern(orig->pattern);
            }
        }

        size_t assignOrAdd = 1; // add
        prepareLoops(handler, hess, evaluations1it, dependentIndexes, assignOrAdd);
    }

}

#endif
