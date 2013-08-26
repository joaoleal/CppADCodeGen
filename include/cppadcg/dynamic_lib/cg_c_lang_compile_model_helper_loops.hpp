#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_INCLUDED
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
     *  Utility classes
     **************************************************************************/

    template<class Base>
    class IndexedDependentLoopInfo {
    public:
        std::vector<size_t> indexes;
        std::vector<CG<Base> > origVals;
        IndexPattern* pattern;

        inline IndexedDependentLoopInfo() :
            pattern(NULL) {
        }
    };

    template<class Base>
    class JacTapeElementLoopInfo {
    public:
        std::map<size_t, IndexedDependentLoopInfo<Base>* > origIndep2Info;
    };

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
        map<LoopAtomicFun<Base>*, map<OperationNode<Base>*, vector<OperationNode<Base>*> > > evaluations1it;

        for (size_t i = 0; i < y.size(); i++) {
            OperationNode<Base>* node = y[i].getOperationNode();
            if (node == NULL || node->getOperationType() != CGLoopAtomicResultOp) {
                continue;
            }

            OperationNode<Base>* loopEvalNode = node->getArguments()[0].getOperation();

#ifndef NDEBUG
            CGOpCode loopOpType = loopEvalNode->getOperationType();
            size_t p = loopEvalNode->getInfo()[3];
            assert(loopOpType == CGLoopForwardOp && p == 0);
#endif

            size_t loopId = loopEvalNode->getInfo()[0];
            LoopAtomicFun<Base>* loop = handler.getLoop(loopId);
            assert(loop != NULL);

            const LoopPosition& depPos = loop->getTapeDependentIndex(i);
            size_t mTape = loop->getTapeDependentCount();
            size_t iteration = depPos.atomic / mTape;

            if (iteration == 0) {
                vector<OperationNode<Base>*>& allLoopEvalResults = evaluations1it[loop][loopEvalNode];
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
        map<LoopAtomicFun<Base>*, map<size_t, map<size_t, JacTapeElementLoopInfo<Base> > > > jacIndexPatterns;

        map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > > dependentIndexes;
        // loop -> loop atomic evaluation -> results
        map<LoopAtomicFun<Base>*, map<OperationNode<Base>*, vector<OperationNode<Base>*> > > evaluations1it;


        for (size_t e = 0; e < nnz; e++) {
            size_t i = _jacSparsity.rows[e];
            size_t j = _jacSparsity.cols[e];
            CGBase& jacVal = jac[e];

            // find LOOP + get loop results
            LoopAtomicFun<Base>* loop = NULL;

            // loop evaluation -> results
            map<OperationNode<Base>*, map<size_t, OperationNode<Base>*> > evals;
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
                OperationNode<Base>* loopEvalNode = evals.begin()->first;

#ifndef NDEBUG
                CGOpCode loopOpType = loopEvalNode->getOperationType();
                size_t p = loopEvalNode->getInfo()[3];
                assert((loopOpType == CGLoopForwardOp && p == 1) ||
                       (loopOpType == CGLoopReverseOp && p == 0));
#endif

                size_t loopId = loopEvalNode->getInfo()[0];
                loop = handler.getLoop(loopId);
                assert(loop != NULL);
            }

            map<size_t, map<size_t, JacTapeElementLoopInfo<Base> > >& lRawIndexes = jacIndexPatterns[loop];

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

            IndexedDependentLoopInfo<Base>* origJacEl = NULL;
            JacTapeElementLoopInfo<Base>& ref = lRawIndexes[tapeI][*tapeJs.begin()];
            if (ref.origIndep2Info.size() != 0) {
                typename map<size_t, IndexedDependentLoopInfo<Base>* >::const_iterator it;
                it = ref.origIndep2Info.find(jRef);
                if (it != ref.origIndep2Info.end()) {
                    origJacEl = it->second;
                }
            }

            if (origJacEl == NULL) {
                // the vector will never allocate more space so this is safe:
                garbageCollection.resize(garbageCollection.size() + 1);
                origJacEl = &garbageCollection.back();
                origJacEl->indexes.resize(loop->getIterationCount(), nnz);
                origJacEl->origVals.resize(loop->getIterationCount());
                ref.origIndep2Info[jRef] = origJacEl;
                dependentIndexes[loop].push_back(origJacEl);
            }
            origJacEl->indexes[iteration] = e;
            origJacEl->origVals[iteration] = jacVal;

            if (iteration == 0) {
                if (evals.size() > 0) {
                    typename map<OperationNode<Base>*, map<size_t, OperationNode<Base>*> >::const_iterator itE;
                    for (itE = evals.begin(); itE != evals.end(); ++itE) {
                        OperationNode<Base>* loopEvalNode = itE->first;
                        const map<size_t, OperationNode<Base>*>& results = itE->second;

                        vector<OperationNode<Base>*>& allLoopEvalResults = evaluations1it[loop][loopEvalNode];

                        size_t result_size = loopEvalNode->getInfo()[5]; // tape result size
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
        typename map<LoopAtomicFun<Base>*, map<size_t, map<size_t, JacTapeElementLoopInfo<Base> > > >::iterator itl;
        for (itl = jacIndexPatterns.begin(); itl != jacIndexPatterns.end(); ++itl) {

            typename map<size_t, map<size_t, JacTapeElementLoopInfo<Base> > >::iterator itI;
            for (itI = itl->second.begin(); itI != itl->second.end(); ++itI) {

                typename map<size_t, JacTapeElementLoopInfo<Base> >::iterator itJ;
                for (itJ = itI->second.begin(); itJ != itI->second.end(); ++itJ) {
                    JacTapeElementLoopInfo<Base>& jacEleInfo = itJ->second;

                    typename map<size_t, IndexedDependentLoopInfo<Base>* >::iterator itO;
                    for (itO = jacEleInfo.origIndep2Info.begin(); itO != jacEleInfo.origIndep2Info.end(); ++itO) {
                        IndexedDependentLoopInfo<Base>* orig = itO->second;

                        // make sure all elements are requested
                        std::vector<size_t>::const_iterator ite;
                        for (ite = orig->indexes.begin(); ite != orig->indexes.end(); ++ite) {
                            if (*ite == nnz) {
                                throw CGException("All jacobian elements of an equation pattern (equation in a loop) must be requested for all iterations");
                            }
                        }

                        orig->pattern = IndexPattern::detect(orig->indexes);
                        handler.manageLoopDependentIndexPattern(orig->pattern);
                    }
                }
            }
        }

        prepareLoops(handler, jac, evaluations1it, dependentIndexes);
    }

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

        /**
         * loops
         */
        size_t mLoopTotal = 0;
        std::vector<bool> eqLoop(m, false);

        std::map<LoopAtomicFun<Base>*, std::vector<CGBase> > loopHess; // hessian elements only for the equation in a given loop

        typename std::set<LoopAtomicFun<Base>* >::const_iterator itl;
        for (itl = _loopAtomics.begin(); itl != _loopAtomics.end(); ++itl) {
            LoopAtomicFun<Base>* loop = *itl;
            size_t iterations = loop->getIterationCount();

            mLoopTotal += loop->getLoopDependentCount();

            vector<CGBase> wLoop(m);
            for (size_t i = 0; i < m; i++) {
                wLoop[i] = Base(0);
            }

            const std::vector<std::vector<LoopPosition> >& depIndexes = loop->getDependentIndexes();
            for (size_t i = 0; i < depIndexes.size(); i++) {
                for (size_t iter = 0; iter < iterations; iter++) {
                    const LoopPosition& pos = depIndexes[i][iter];
                    eqLoop[pos.original] = true;
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

        if (mLoopTotal < m) {
            /**
             * equations not in loops
             */
            vector<CGBase> ww(m);
            for (size_t i = 0; i < m; i++) {
                ww[i] = eqLoop[i] ? Base(0) : w[i];
            }

            CppAD::sparse_hessian_work work;
            vector<CGBase> lowerHess(lowerHessRows.size());
            _fun->SparseHessian(indVars, ww, _hessSparsity.sparsity, lowerHessRows, lowerHessCols, lowerHess, work);

            std::vector<CGBase> hessNoLoops(_hessSparsity.rows.size());
            for (size_t i = 0; i < lowerHessOrder.size(); i++) {
                hessNoLoops[lowerHessOrder[i]] = lowerHess[i];
            }

            // make use of the symmetry of the Hessian in order to reduce operations
            std::map<size_t, size_t>::const_iterator it2;
            for (it2 = duplicates.begin(); it2 != duplicates.end(); ++it2) {
                hessNoLoops[it2->first] = hessNoLoops[it2->second];
            }

            for (size_t e = 0; e < hess.size(); e++) {
                hess[e] += hessNoLoops[e];
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
        map<size_t, map<size_t, map<size_t, map<size_t, IndexedDependentLoopInfo<Base>* > > > > hessIndexPatterns;

        map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > > dependentIndexes;
        // loop -> loop atomic evaluation -> results
        map<LoopAtomicFun<Base>*, map<OperationNode<Base>*, vector<OperationNode<Base>*> > > evaluations1it;

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
            map<OperationNode<Base>*, map<size_t, OperationNode<Base>*> > evals;
            findLoopEvaluations(handler, loop, hessVal.getOperationNode(), evals); ////////////

            if (evals.empty())
                continue;

            if (evals.size() > 1 || evals.begin()->second.size() > 1) {
                throw CGException("Unable to generate expression for an hessian element which "
                                  "is either associated with multiple indexed independents "
                                  "or an indepedent with constant index.");
            }
            /*
                        typename map<OperationNode<Base>*, map<size_t, OperationNode<Base>*> >::const_iterator itEval;
                        for (itEval = evals.begin(); itEval != evals.end(); ++itEval) {
                            OperationNode<Base>* loopEvalNode = itEval->first;
                            CGOpCode loopOpType = loopEvalNode->getOperationType();
                            size_t p = loopEvalNode->getInfo()[3];

                            if (loopOpType == CGLoopReverseOp && p == 1) {
                                typename map<size_t, OperationNode<Base>*>::const_iterator itRes;
                                for (itRes = itEval->second.begin(); itRes != itEval->second.end(); ++itRes) {
                                    size_t jAtomic = itRes->first / 2;

                                    const std::map<size_t, size_t>& tapeJ2iter = loop->getAtomicIndependentLocations(jAtomic);
                                    std::map<size_t, size_t>::const_iterator t2i;
                                    for (t2i = tapeJ2iter.begin(); t2i != tapeJ2iter.end(); ++t2i) {
                                        // jTape can appear only linearly in the model and 
                                        // therefore it might not influence the hessian element
                                        size_t jTape = t2i->first;
                                        size_t iteration = t2i->second;

                                    }

                                    if (!loop->isIndexed(jTape)) {
                                        throw CGException("Found an hessian element which is not indexed! "
                                                          "Hessian elements must correspond to indexed independent variables.");
                                    }


                                    //iteration = loop->getIterationOfIndexedIndep(jAtomic, origJ);
                                }
                            }
                        }
             */
#ifndef NDEBUG
            {
                OperationNode<Base>* loopEvalNode = evals.begin()->first;
                CGOpCode loopOpType = loopEvalNode->getOperationType();
                size_t loopId = loopEvalNode->getInfo()[0];
                size_t p = loopEvalNode->getInfo()[3];
                assert(loopOpType == CGLoopReverseOp && p == 1);
                assert(loop == handler.getLoop(loopId));
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


            set<size_t> tapeJ2s = loop->getIndependentTapeIndexes(j2);
            assert(tapeJ2s.size() >= 1); // if >1 then the independent is used by several temporary variables

            bool isTemporary2 = loop->isTemporary(*tapeJ2s.begin());
            size_t jRef2 = isTemporary2 ? j2 : 0;

            IndexedDependentLoopInfo<Base>* origHessEl = NULL;
            map<size_t, IndexedDependentLoopInfo<Base>* >& ref = hessIndexPatterns[*tapeJ1s.begin()][jRef1][*tapeJ2s.begin()];
            if (ref.size() != 0) {
                typename map<size_t, IndexedDependentLoopInfo<Base>* >::const_iterator it;
                it = ref.find(jRef2);
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
                ref[jRef2] = origHessEl;
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
                typename map<OperationNode<Base>*, map<size_t, OperationNode<Base>*> >::const_iterator itE;
                for (itE = evals.begin(); itE != evals.end(); ++itE) {
                    OperationNode<Base>* loopEvalNode = itE->first;
                    const map<size_t, OperationNode<Base>*>& results = itE->second;

                    size_t result_size = loopEvalNode->getInfo()[5]; // tape result size
#ifndef NDEBUG
                    size_t loopId = loopEvalNode->getInfo()[0];
                    assert(handler.getLoop(loopId) == loop);
#endif
                    vector<OperationNode<Base>*>& allLoopEvalResults = evaluations1it[loop][loopEvalNode];
                    allLoopEvalResults.resize(result_size);

                    typename map<size_t, OperationNode<Base>*>::const_iterator itR;
                    for (itR = results.begin(); itR != results.end(); ++itR) {
                        size_t tapeIndex = itR->second->getInfo()[1];
                        allLoopEvalResults[tapeIndex] = itR->second;
                    }
                }
            }
        }

        /**
         * Generate index patterns for the dependent variables
         */
        typename map<size_t, map<size_t, map<size_t, map<size_t, IndexedDependentLoopInfo<Base>* > > > >::iterator itJ1;
        for (itJ1 = hessIndexPatterns.begin(); itJ1 != hessIndexPatterns.end(); ++itJ1) {

            typename map<size_t, map<size_t, map<size_t, IndexedDependentLoopInfo<Base>* > > >::iterator itO1;
            for (itO1 = itJ1->second.begin(); itO1 != itJ1->second.end(); ++itO1) {

                typename map<size_t, map<size_t, IndexedDependentLoopInfo<Base>* > >::iterator itJ2;
                for (itJ2 = itO1->second.begin(); itJ2 != itO1->second.end(); ++itJ2) {

                    typename map<size_t, IndexedDependentLoopInfo<Base>* >::iterator itO2;
                    for (itO2 = itJ2->second.begin(); itO2 != itJ2->second.end(); ++itO2) {
                        IndexedDependentLoopInfo<Base>* orig = itO2->second;

                        // make sure all element are requested
                        std::vector<size_t>::const_iterator ite;
                        for (ite = orig->indexes.begin(); ite != orig->indexes.end(); ++ite) {
                            if (*ite == nnz) {
                                throw CGException("All hessian elements must be requested for all iterations");
                            }
                        }

                        orig->pattern = IndexPattern::detect(orig->indexes);
                        handler.manageLoopDependentIndexPattern(orig->pattern);
                    }
                }
            }
        }

        size_t assignOrAdd = 1; // add
        prepareLoops(handler, hess, evaluations1it, dependentIndexes, assignOrAdd);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareLoops(CodeHandler<Base>& handler,
                                                     std::vector<CGBase>& dep,
                                                     std::map<LoopAtomicFun<Base>*, std::map<OperationNode<Base>*, vector<OperationNode<Base>*> > >& evaluations1it,
                                                     std::map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > >& dependentIndexes,
                                                     size_t assignOrAdd) {

        std::map<LoopAtomicFun<Base>*, vector<OperationNode<Base>* > > loopIndexedIndependents;
        std::map<LoopAtomicFun<Base>*, vector<OperationNode<Base>* > > loopIndexedIndependents2;

        /**
         * insert the loops into the operation graph
         */
        typename std::map<LoopAtomicFun<Base>*, std::map<OperationNode<Base>*, vector<OperationNode<Base>*> > >::iterator itL;
        for (itL = evaluations1it.begin(); itL != evaluations1it.end(); ++itL) {
            LoopAtomicFun<Base>* loopFunc = itL->first;
            std::map<OperationNode<Base>*, vector<OperationNode<Base>*> >& loopAtomEvaluations1it = itL->second;

            std::vector<size_t> startEndInfo(2);
            startEndInfo[0] = loopFunc->getLoopId();
            startEndInfo[1] = loopFunc->getIterationCount();
            std::vector<Argument<Base> > startArgs;
            OperationNode<Base>* loopStart = new OperationNode<Base>(CGLoopStartOp, startEndInfo, startArgs);
            handler.manageOperationNodeMemory(loopStart);

            // indexed independents
            vector<OperationNode<Base>* >& indexedIndependents = loopIndexedIndependents[loopFunc]; // zero order
            vector<OperationNode<Base>* >& indexedIndependents2 = loopIndexedIndependents2[loopFunc]; // first order
            if (indexedIndependents.size() == 0) {
                size_t nIndexed = loopFunc->getIndependentIndexPatterns().size();
                indexedIndependents.resize(nIndexed); // zero order
                std::vector<Argument<Base> > xIndexedArgs(1);
                xIndexedArgs[0] = Argument<Base>(*loopStart);
                std::vector<size_t> info(2);
                info[0] = 0; // tx
                for (size_t j = 0; j < nIndexed; j++) {
                    info[1] = j;
                    indexedIndependents[j] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, xIndexedArgs);
                    handler.manageOperationNodeMemory(indexedIndependents[j]);
                }
            }

            std::vector<Argument<Base> > endArgs;

            typename std::map<OperationNode<Base>*, vector<OperationNode<Base>*> >::iterator itE;
            for (itE = loopAtomEvaluations1it.begin(); itE != loopAtomEvaluations1it.end(); ++itE) {
                OperationNode<Base>* loopEvalNode = itE->first;
                vector<OperationNode<Base>*>& resultNodes1it = itE->second;

                // evaluate tape 1st iteration (generate the operation graph)
                if (loopEvalNode != NULL) {
                    vector<CG<Base> > newTapeResults1it = evalLoopTape(handler, *loopFunc, *loopEvalNode,
                                                                       indexedIndependents, indexedIndependents2,
                                                                       loopStart);

                    assert(resultNodes1it.size() == newTapeResults1it.size());
                    for (size_t j = 0; j < resultNodes1it.size(); j++) {
                        if (resultNodes1it[j] != NULL) {
                            assert(resultNodes1it[j]->getOperationType() == CGLoopAtomicResultOp);
                            resultNodes1it[j]->getArguments()[0] = asArgument(newTapeResults1it[j]);
                        }
                    }
                }
            }

            /**
             * change the dependents
             */
            const vector<IndexedDependentLoopInfo<Base>* >& dependents = dependentIndexes.at(loopFunc);
            size_t dep_size = dependents.size();
            std::vector<Argument<Base> > indexedArgs(1);
            std::vector<size_t> info(2);
            for (size_t i = 0; i < dep_size; i++) {
                IndexedDependentLoopInfo<Base>& depInfo = *dependents[i];

                indexedArgs[0] = asArgument(depInfo.origVals[0]); // value from first iteration!
                info[0] = handler.addLoopDependentIndexPattern(*depInfo.pattern); // dependent index pattern location
                info[1] = assignOrAdd;

                OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, info, indexedArgs);
                handler.manageOperationNodeMemory(yIndexed);
                endArgs.push_back(Argument<Base>(*yIndexed));
            }

            OperationNode<Base>* loopEnd = new OperationNode<Base>(CGLoopEndOp, startEndInfo, endArgs);
            handler.manageOperationNodeMemory(loopEnd);

            for (size_t i = 0; i < dep_size; i++) {
                for (size_t it = 0; it < dependents[i]->indexes.size(); it++) {
                    size_t e = dependents[i]->indexes[it];
                    dep[e] = handler.createCG(Argument<Base>(*loopEnd));
                }
            }

            moveNonIndexedOutsideLoop(*loopStart, *loopEnd);
        }

    }

    template<class Base>
    void CLangCompileModelHelper<Base>::findLoopEvaluations(CodeHandler<Base>& handler,
                                                            OperationNode<Base>* node,
                                                            std::map<OperationNode<Base>*, std::map<size_t, OperationNode<Base>*> >& evals) {
        if (node == NULL) {
            return;
        }

        CGOpCode op = node->getOperationType();
        if (op == CGLoopAtomicResultOp) {
            size_t pos = node->getInfo()[0];
            OperationNode<Base>* loopEvalNode = node->getArguments()[0].getOperation();
            evals[loopEvalNode][pos] = node;
        } else {
            const std::vector<Argument<Base> >& args = node->getArguments();
            for (size_t a = 0; a < args.size(); a++) {
                findLoopEvaluations(handler, args[a].getOperation(), evals);
            }
        }

    }

    template<class Base>
    void CLangCompileModelHelper<Base>::findLoopEvaluations(CodeHandler<Base>& handler,
                                                            LoopAtomicFun<Base>* loop,
                                                            OperationNode<Base>* node,
                                                            std::map<OperationNode<Base>*, std::map<size_t, OperationNode<Base>*> >& evals) {
        if (node == NULL) {
            return;
        }

        CGOpCode op = node->getOperationType();
        if (op == CGLoopAtomicResultOp) {
            size_t pos = node->getInfo()[0];
            OperationNode<Base>* loopEvalNode = node->getArguments()[0].getOperation();
            size_t loopId = loopEvalNode->getInfo()[0];
            if (loopId == loop->getLoopId()) {
                evals[loopEvalNode][pos] = node;
            }
        } else {
            const std::vector<Argument<Base> >& args = node->getArguments();
            for (size_t a = 0; a < args.size(); a++) {
                findLoopEvaluations(handler, loop, args[a].getOperation(), evals);
            }
        }
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::evalLoopTape(CodeHandler<Base>& handler,
                                                                  LoopAtomicFun<Base>& atomic,
                                                                  const OperationNode<Base>& loopEvalNode,
                                                                  const vector<OperationNode<Base>*>& indexedIndependents,
                                                                  vector<OperationNode<Base>* >& indexedIndependents2,
                                                                  OperationNode<Base>* loopStart) {
        assert(loopEvalNode.getOperationType() == CGLoopForwardOp || loopEvalNode.getOperationType() == CGLoopReverseOp);
        assert(loopEvalNode.getInfo().size() == 6);

        size_t p = loopEvalNode.getInfo()[3];

        if (loopEvalNode.getOperationType() == CGLoopForwardOp) {
            /**
             * Forward mode
             */
            if (p == 0) {
                return generateLoopForward0Graph(handler, atomic, indexedIndependents, loopEvalNode.getArguments());
            } else if (p == 1) {
                return generateForward1Graph(handler, atomic,
                                             indexedIndependents, indexedIndependents2,
                                             loopEvalNode.getArguments(), loopStart);
            }
        } else {
            /**
             * Reverse mode
             */
            if (p == 0) {
                return generateReverse1Graph(handler, atomic,
                                             indexedIndependents,
                                             loopEvalNode.getArguments(), loopStart);
            } else if (p == 1) {
                return generateReverse2Graph(handler, atomic,
                                             indexedIndependents, indexedIndependents2,
                                             loopEvalNode.getArguments(), loopStart);
            }
        }

        throw CGException("Unable to generate operation graph for loop: not implemented yet!");
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::generateLoopForward0Graph(CodeHandler<Base>& handler,
                                                                               LoopAtomicFun<Base>& atomic,
                                                                               const vector<OperationNode<Base>*>& indexedIndependents,
                                                                               const std::vector<Argument<Base> >& args) {
        typedef CppAD::CG<Base> CGB;

        ADFun<CGB>* fun = atomic.getTape();

        vector<CGB> tx = createLoopGraphIndependentVector(handler, atomic, indexedIndependents, args, 0);
        vector<CGB> ty = fun->Forward(0, tx);

        return ty;
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::generateForward1Graph(CodeHandler<Base>& handler,
                                                                           LoopAtomicFun<Base>& atomic,
                                                                           const vector<OperationNode<Base>*>& indexedIndependents,
                                                                           vector<OperationNode<Base>*>& indexedIndependents2,
                                                                           const std::vector<Argument<Base> >& argsAtomic, //argsAtomic = [txAtomic]
                                                                           OperationNode<Base>* loopStart) {
        typedef CppAD::CG<Base> CGB;

        ADFun<CGB>* fun = atomic.getTape();

        // zero order
        vector<CGB> x = createLoopGraphIndependentVector(handler, atomic, indexedIndependents, argsAtomic, 1);
        assert(x.size() == fun->Domain());

        fun->Forward(0, x);

        // forward first order
        vector<CGB> tx = createLoopGraphIndependentVectorTx2(handler, atomic, x, indexedIndependents2, argsAtomic, 1, loopStart);
        assert(tx.size() == fun->Domain() * 2);

        vector<CGB> ty = fun->Forward(1, tx);

        return ty;
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::generateReverse1Graph(CodeHandler<Base>& handler,
                                                                           LoopAtomicFun<Base>& atomic,
                                                                           const vector<OperationNode<Base>*>& indexedIndependents,
                                                                           const std::vector<Argument<Base> >& argsAtomic, //argsAtomic = [txAtomic pyAtomic]
                                                                           OperationNode<Base>* loopStart) {
        typedef CppAD::CG<Base> CGB;

        size_t m = atomic.getTapeDependentCount();
        size_t nFull = atomic.getLoopIndependentCount();

        ADFun<CGB>* fun = atomic.getTape();

        // zero order
        vector<CGB> x = createLoopGraphIndependentVector(handler, atomic, indexedIndependents, argsAtomic, 0);
        assert(x.size() == fun->Domain());

        fun->Forward(0, x);

        // reverse first order
        vector<CGB> py(m);
        for (size_t i = 0; i < m; i++) {
            if (argsAtomic[nFull + i].getOperation() == NULL) { // first iteration
                py[i] = handler.createCG(argsAtomic[nFull + i]);
            } else {
                std::vector<size_t> info(2);
                info[0] = 1; // py
                info[1] = i;
                std::vector<Argument<Base> > args(1);
                args[0] = *loopStart;
                OperationNode<Base>* indexedPy = new OperationNode<Base>(CGLoopIndexedIndepOp, info, args);
                py[i] = handler.createCG(indexedPy);
            }
        }

        vector<CGB> px = fun->Reverse(1, py);

        return px;
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::generateReverse2Graph(CodeHandler<Base>& handler,
                                                                           LoopAtomicFun<Base>& atomic,
                                                                           const vector<OperationNode<Base>*>& indexedIndependents,
                                                                           vector<OperationNode<Base>*>& indexedIndependents2,
                                                                           const std::vector<Argument<Base> >& argsAtomic, //argsAtomic = [txAtomic pyAtomic]
                                                                           OperationNode<Base>* loopStart) {
        typedef CppAD::CG<Base> CGB;

        size_t m = atomic.getTapeDependentCount();
        size_t nFull = atomic.getLoopIndependentCount();

        ADFun<CGB>* fun = atomic.getTape();

        // zero order
        vector<CGB> x = createLoopGraphIndependentVector(handler, atomic, indexedIndependents, argsAtomic, 1);
        assert(x.size() == fun->Domain());

        fun->Forward(0, x);

        // forward first order
        vector<CGB> tx = createLoopGraphIndependentVectorTx2(handler, atomic, x, indexedIndependents2, argsAtomic, 1, loopStart);
        assert(tx.size() == fun->Domain() * 2);

        fun->Forward(1, tx);

        // reverse second order
        vector<CGB> py(2 * m);
        for (size_t i = 0; i < 2 * m; i++) {
            if (argsAtomic[nFull * 2 + i].getOperation() == NULL) {
                py[i] = handler.createCG(argsAtomic[nFull * 2 + i]);
            } else {
                std::vector<size_t> info(2);
                info[0] = 1; // py
                info[1] = i / 2; // TODO: improve this
                std::vector<Argument<Base> > args(1);
                args[0] = *loopStart;
                OperationNode<Base>* indexedPy = new OperationNode<Base>(CGLoopIndexedIndepOp, info, args);
                py[i] = handler.createCG(indexedPy);
            }
        }

        vector<CGB> px = fun->Reverse(2, py);

        return px;
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::createLoopGraphIndependentVector(CodeHandler<Base>& handler,
                                                                                      LoopAtomicFun<Base>& atomic,
                                                                                      const vector<OperationNode<Base>*>& indexedIndependents,
                                                                                      const std::vector<Argument<Base> >& argsAtomic,
                                                                                      size_t p) {
        typedef CppAD::CG<Base> CGB;

        const std::vector<LoopPositionTmp>& temporaryIndependents = atomic.getTemporaryIndependents();
        const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = atomic.getIndexedIndepIndexes();
        const std::vector<LoopPosition>& nonIndexedIndepIndexes = atomic.getNonIndexedIndepIndexes();

        size_t nIndexed = indexedIndepIndexes.size();
        size_t nNonIndexed = nonIndexedIndepIndexes.size();
        size_t nTape = nIndexed + nNonIndexed + temporaryIndependents.size();

        // indexed independents
        vector<CGB> x(nTape);
        for (size_t j = 0; j < nIndexed; j++) {
            x[j] = handler.createCG(Argument<Base>(*indexedIndependents[j]));
        }
        // non indexed
        for (size_t j = 0; j < nNonIndexed; j++) {
            x[nIndexed + j] = handler.createCG(argsAtomic[nonIndexedIndepIndexes[j].atomic * (p + 1)]);
        }
        // temporaries
        for (size_t j = 0; j < temporaryIndependents.size(); j++) {
            x[nIndexed + nNonIndexed + j] = handler.createCG(argsAtomic[temporaryIndependents[j].atomic * (p + 1)]);
        }
        return x;
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::createLoopGraphIndependentVectorTx2(CodeHandler<Base>& handler,
                                                                                         LoopAtomicFun<Base>& atomic,
                                                                                         const vector<CG<Base> >& x,
                                                                                         vector<OperationNode<Base>*>& indexedIndependents2,
                                                                                         const std::vector<Argument<Base> >& argsAtomic,
                                                                                         size_t p,
                                                                                         OperationNode<Base>* loopStart) {
        typedef CppAD::CG<Base> CGB;

        const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = atomic.getIndexedIndepIndexes();
        const std::vector<LoopPosition>& nonIndexedIndepIndexes = atomic.getNonIndexedIndepIndexes();
        const std::vector<LoopPositionTmp>& temporaryIndependents = atomic.getTemporaryIndependents();

        if (indexedIndependents2.size() == 0) {
            size_t nIndexed = indexedIndepIndexes.size();
            indexedIndependents2.resize(nIndexed); // first order
            std::vector<Argument<Base> > xIndexedArgs(1);
            xIndexedArgs[0] = Argument<Base>(*loopStart);
            std::vector<size_t> info(2);
            info[0] = 1; // tx1
            for (size_t j = 0; j < nIndexed; j++) {
                info[1] = j;
                indexedIndependents2[j] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, xIndexedArgs);
                handler.manageOperationNodeMemory(indexedIndependents2[j]);
            }
        }

        ADFun<CGB>* fun = atomic.getTape();
        size_t nTape = fun->Domain();
        size_t nIndexed = indexedIndepIndexes.size();
        size_t nNonIndexed = nonIndexedIndepIndexes.size();
        size_t nTmp = temporaryIndependents.size();

        vector<CGB> tx(nTape * 2);
        for (size_t j = 0; j < nTape; j++) {
            tx[j * 2] = x[j];
        }

        size_t k1 = p + 1;

        for (size_t j = 0; j < nIndexed; j++) {
            const LoopPosition& pos = indexedIndepIndexes[j][0];
            if (argsAtomic[pos.atomic * k1 + 1].getOperation() == NULL) {
                tx[pos.tape * k1 + 1] = handler.createCG(argsAtomic[pos.atomic * k1 + 1]);
            } else {
                tx[pos.tape * k1 + 1] = handler.createCG(Argument<Base>(*indexedIndependents2[j]));
            }
        }
        for (size_t j = 0; j < nNonIndexed; j++) {
            const LoopPosition& pos = nonIndexedIndepIndexes[j];
            tx[pos.tape * k1 + 1] = handler.createCG(argsAtomic[pos.atomic * k1 + 1]);
        }
        for (size_t j = 0; j < nTmp; j++) {
            const LoopPositionTmp& pos = temporaryIndependents[j];
            tx[pos.tape * k1 + 1] = handler.createCG(argsAtomic[pos.atomic * k1 + 1]);
        }

        return tx;
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::moveNonIndexedOutsideLoop(OperationNode<Base>& loopStart,
                                                                  OperationNode<Base>& loopEnd) {
        //EquationPattern<Base>::uncolor(dependents[dep].getOperationNode());
        std::set<OperationNode<Base>*> nonIndexed;

        const std::vector<Argument<Base> >& endArgs = loopEnd.getArguments();
        for (size_t i = 0; i < endArgs.size(); i++) {
            assert(endArgs[i].getOperation() != NULL);
            findNonIndexedNodes(*endArgs[i].getOperation(), nonIndexed);
        }

        std::vector<Argument<Base> >& startArgs = loopStart.getArguments();
        assert(startArgs.empty());

        startArgs.resize(nonIndexed.size());
        size_t i = 0;
        typename std::set<OperationNode<Base>*>::const_iterator it;
        for (it = nonIndexed.begin(); it != nonIndexed.end(); ++it, i++) {
            startArgs[i] = Argument<Base>(**it);
        }
    }

    template<class Base>
    bool CLangCompileModelHelper<Base>::findNonIndexedNodes(OperationNode<Base>& node,
                                                            std::set<OperationNode<Base>*>& nonIndexed) {
        if (node.getColor() > 0)
            return node.getColor() == 1;

        if (node.getOperationType() == CGLoopIndexedIndepOp) {
            node.setColor(2);
            return false; // does NOT depend on an index
        }

        const std::vector<Argument<Base> >& args = node.getArguments();
        size_t size = args.size();

        bool indexedPath = false; // whether or not this node depends on indexed independents
        bool nonIndexedArgs = false; // whether or not there are non indexed arguments
        for (size_t a = 0; a < size; a++) {
            OperationNode<Base>* arg = args[a].getOperation();
            if (arg != NULL) {
                bool nonIndexedArg = findNonIndexedNodes(*arg, nonIndexed);
                nonIndexedArgs |= nonIndexedArg;
                indexedPath |= !nonIndexedArg;
            }
        }

        node.setColor(indexedPath ? 2 : 1);

        if (indexedPath && nonIndexedArgs) {
            for (size_t a = 0; a < size; a++) {
                OperationNode<Base>* arg = args[a].getOperation();
                if (arg != NULL && arg->getColor() == 1 && arg->getOperationType() != CGInvOp) {
                    nonIndexed.insert(arg);
                }
            }
        }

        return !indexedPath;
    }
}

#endif
