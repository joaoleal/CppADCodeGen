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
        std::vector<size_t> jacIndexes;
        std::vector<CG<Base> > origJacVals;
        IndexPattern* jacPattern;

        inline IndexedDependentLoopInfo() :
            jacPattern(NULL) {
        }
    };

    template<class Base>
    class JacOrigElementIndepLoopInfo {
    public:
        IndexedDependentLoopInfo<Base>* origJacElement; // original model jacobian element

        inline JacOrigElementIndepLoopInfo() :
            origJacElement(NULL) {
        }
    };

    template<class Base>
    class JacTapeElementLoopInfo {
    public:
        std::map<size_t, JacOrigElementIndepLoopInfo<Base> > origIndep2Info;
    };

    /***************************************************************************
     *  Methods related with loop insertion into the operation graph
     **************************************************************************/

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareForward0WithLoops(CodeHandler<Base>& handler,
                                                                 std::vector<CGBase>& y) {

        std::vector<IndexedDependentLoopInfo<Base> > garbageCollection(y.size()); ////////// <- rethink this!!!!!!!
        garbageCollection.resize(y.size());

        std::map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > > dependentIndexes;
        // loop -> loop atomic evaluation -> results
        std::map<LoopAtomicFun<Base>*, std::map<OperationNode<Base>*, vector<OperationNode<Base>*> > > evaluations1it;

        for (size_t i = 0; i < y.size(); i++) {
            OperationNode<Base>* node = y[i].getOperationNode();
            if (node == NULL || node->getOperationType() != CGLoopAtomicResultOp) {
                continue;
            }

            OperationNode<Base>* loopEvalNode = node->getArguments()[0].getOperation();

            CGOpCode loopOpType = loopEvalNode->getOperationType();
            size_t loopId = loopEvalNode->getInfo()[0];
            //size_t p = loopNode->getInfo()[3];
            assert(loopOpType == CGLoopForwardOp && loopEvalNode->getInfo()[3] == 0);

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
            origEl->jacIndexes.resize(loop->getIterationCount());
            origEl->origJacVals.resize(loop->getIterationCount());
            origEl->jacIndexes[iteration] = i;
            origEl->origJacVals[iteration] = y[i];
            origEl->jacPattern = loop->getDependentIndexPatterns()[depPos.tape];
        }

        prepareLoops(handler, y, evaluations1it, dependentIndexes);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareSparseJacobianWithLoops(CodeHandler<Base>& handler,
                                                                       std::vector<CGBase>& jac) {
        size_t nnz = _jacSparsity.rows.size();

        std::vector<IndexedDependentLoopInfo<Base> > garbageCollection; ////////// <- rethink this!!!!!!!
        garbageCollection.reserve(nnz);

        /**
         * Generate index patterns for the jacobian elements resulting from loops
         */
        // loop -> equation pattern -> tape independent -> orig independent(temporaries only) -> iteration = position
        std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > > jacIndexPatterns;

        std::map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > > dependentIndexes;
        // loop -> loop atomic evaluation -> results
        std::map<LoopAtomicFun<Base>*, std::map<OperationNode<Base>*, vector<OperationNode<Base>*> > > evaluations1it;


        for (size_t e = 0; e < nnz; e++) {
            // find LOOP + get loop results

            // loop evaluation -> results
            std::map<OperationNode<Base>*, std::map<size_t, OperationNode<Base>*> > evals;
            findLoopEvaluations(handler, jac[e].getOperationNode(), evals);

            if (evals.empty())
                continue;

            assert(evals.size() == 1);
            OperationNode<Base>* loopEvalNode = evals.begin()->first;

            CGOpCode loopOpType = loopEvalNode->getOperationType();
            size_t loopId = loopEvalNode->getInfo()[0];
            //size_t p = loopNode->getInfo()[3];
            assert((loopOpType == CGLoopForwardOp && loopEvalNode->getInfo()[3] == 1) ||
                   (loopOpType == CGLoopReverseOp && loopEvalNode->getInfo()[3] == 0));

            LoopAtomicFun<Base>* loop = handler.getLoop(loopId);
            assert(loop != NULL);
            std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > >& lRawIndexes = jacIndexPatterns[loop];

            size_t i = _jacSparsity.rows[e];
            size_t j = _jacSparsity.cols[e];

            const LoopPosition& depPos = loop->getTapeDependentIndex(i);
            size_t tapeI = depPos.tape;
            if (!loop->isIndependentVariablesFullyDetached(tapeI)) {
                throw CGException("There are independent variables which appear as indexed and non-indexed for the same equation pattern");
            }
            size_t iteration = depPos.atomic / loop->getTapeDependentCount();

            // must get a single tape independent index!!!!!
            std::set<size_t> tapeJs = loop->getIndependentTapeIndexes(tapeI, j);
            assert(tapeJs.size() >= 1); // if >1 then the independent is used by several temporary variables

            bool isTemporary = loop->isTemporary(*tapeJs.begin());
            size_t jRef = isTemporary ? j : 0;

            IndexedDependentLoopInfo<Base>* origJacEl = NULL;
            JacTapeElementLoopInfo<Base>& ref = lRawIndexes[tapeI][*tapeJs.begin()];
            if (ref.origIndep2Info.size() != 0) {
                typename std::map<size_t, JacOrigElementIndepLoopInfo<Base> >::const_iterator it;
                it = ref.origIndep2Info.find(jRef);
                if (it != ref.origIndep2Info.end()) {
                    origJacEl = it->second.origJacElement;
                }
            }

            if (origJacEl == NULL) {
                // the vector will never allocate more space so this is safe:
                garbageCollection.resize(garbageCollection.size() + 1);
                origJacEl = &garbageCollection.back();
                origJacEl->jacIndexes.resize(loop->getIterationCount(), nnz);
                origJacEl->origJacVals.resize(loop->getIterationCount());
                ref.origIndep2Info[jRef].origJacElement = origJacEl;
                dependentIndexes[loop].push_back(origJacEl);
            }
            origJacEl->jacIndexes[iteration] = e;
            origJacEl->origJacVals[iteration] = jac[e];

            if (iteration == 0) {
                typename std::map<OperationNode<Base>*, std::map<size_t, OperationNode<Base>*> >::const_iterator itE;
                for (itE = evals.begin(); itE != evals.end(); ++itE) {
                    OperationNode<Base>* loopEvalNode = itE->first;
                    const std::map<size_t, OperationNode<Base>*>& results = itE->second;

                    size_t loopId = loopEvalNode->getInfo()[0];
                    size_t result_size = loopEvalNode->getInfo()[5]; // tape result size
                    LoopAtomicFun<Base>* loop = handler.getLoop(loopId);
                    assert(loop != NULL);

                    vector<OperationNode<Base>*>& allLoopEvalResults = evaluations1it[loop][loopEvalNode];
                    allLoopEvalResults.resize(result_size);

                    typename std::map<size_t, OperationNode<Base>*>::const_iterator itR;
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
        typename std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > >::iterator itl;
        for (itl = jacIndexPatterns.begin(); itl != jacIndexPatterns.end(); ++itl) {

            typename std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > >::iterator itI;
            for (itI = itl->second.begin(); itI != itl->second.end(); ++itI) {

                typename std::map<size_t, JacTapeElementLoopInfo<Base> >::iterator itJ;
                for (itJ = itI->second.begin(); itJ != itI->second.end(); ++itJ) {
                    JacTapeElementLoopInfo<Base>& jacEleInfo = itJ->second;

                    typename std::map<size_t, JacOrigElementIndepLoopInfo<Base> >::iterator itO;
                    for (itO = jacEleInfo.origIndep2Info.begin(); itO != jacEleInfo.origIndep2Info.end(); ++itO) {
                        JacOrigElementIndepLoopInfo<Base>& tapeEleInfo = itO->second;
                        IndexedDependentLoopInfo<Base>* orig = tapeEleInfo.origJacElement;

                        // make sure all element are requested
                        std::vector<size_t>::const_iterator ite;
                        for (ite = orig->jacIndexes.begin(); ite != orig->jacIndexes.end(); ++ite) {
                            if (*ite == nnz) {
                                throw CGException("All jacobian elements of an equation pattern (equation in a loop) must be requested for all iterations");
                            }
                        }

                        orig->jacPattern = IndexPattern::detect(orig->jacIndexes);
                        handler.manageLoopDependentIndexPattern(orig->jacPattern);
                    }
                }
            }
        }

        prepareLoops(handler, jac, evaluations1it, dependentIndexes);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareLoops(CodeHandler<Base>& handler,
                                                     std::vector<CGBase>& dep,
                                                     std::map<LoopAtomicFun<Base>*, std::map<OperationNode<Base>*, vector<OperationNode<Base>*> > >& evaluations1it,
                                                     std::map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > >& dependentIndexes) {
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
            std::vector<Argument<Base> > startArgs; ////////////////////// must find non-indexed nodes
            OperationNode<Base>* loopStart = new OperationNode<Base>(CGLoopStartOp, startEndInfo, startArgs);
            handler.manageOperationNodeMemory(loopStart);

            // indexed independents
            size_t nIndexed = loopFunc->getIndependentIndexPatterns().size();
            vector<OperationNode<Base>* > indexedIndependents(nIndexed);
            std::vector<Argument<Base> > xIndexedArgs(1);
            xIndexedArgs[0] = Argument<Base>(*loopStart);
            std::vector<size_t> info(2);
            info[0] = 0;
            for (size_t j = 0; j < nIndexed; j++) {
                info[1] = j;
                indexedIndependents[j] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, xIndexedArgs);
                handler.manageOperationNodeMemory(indexedIndependents[j]);
            }

            std::vector<Argument<Base> > endArgs;

            typename std::map<OperationNode<Base>*, vector<OperationNode<Base>*> >::iterator itE;
            for (itE = loopAtomEvaluations1it.begin(); itE != loopAtomEvaluations1it.end(); ++itE) {
                OperationNode<Base>* loopEvalNode = itE->first;
                vector<OperationNode<Base>*>& resultNodes1it = itE->second;

                // evaluate tape 1st iteration (generate the operation graph)
                vector<CG<Base> > newTapeResults1it = evalLoopTape(handler, *loopFunc, *loopEvalNode, indexedIndependents);

                assert(resultNodes1it.size() == newTapeResults1it.size());
                for (size_t j = 0; j < resultNodes1it.size(); j++) {
                    if (resultNodes1it[j] != NULL) {
                        assert(resultNodes1it[j]->getOperationType() == CGLoopAtomicResultOp);
                        resultNodes1it[j]->getArguments()[0] = asArgument(newTapeResults1it[j]);
                    }
                }
            }

            /**
             * change the dependents
             */
            const vector<IndexedDependentLoopInfo<Base>* >& dependents = dependentIndexes.at(loopFunc);
            size_t dep_size = dependents.size();
            std::vector<Argument<Base> > indexedArgs(1);
            for (size_t i = 0; i < dep_size; i++) {
                IndexedDependentLoopInfo<Base>& depInfo = *dependents[i];

                assert(depInfo.origJacVals[0].getOperationNode() != NULL);
                indexedArgs[0] = Argument<Base>(*depInfo.origJacVals[0].getOperationNode()); // value from first iteration!
                info[0] = handler.addLoopDependentIndexPattern(*depInfo.jacPattern); // dependent index pattern location

                OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, info, indexedArgs);
                handler.manageOperationNodeMemory(yIndexed);
                endArgs.push_back(Argument<Base>(*yIndexed));
            }

            OperationNode<Base>* loopEnd = new OperationNode<Base>(CGLoopEndOp, startEndInfo, endArgs);
            handler.manageOperationNodeMemory(loopEnd);

            for (size_t i = 0; i < dep_size; i++) {
                for (size_t it = 0; it < dependents[i]->jacIndexes.size(); it++) {
                    size_t e = dependents[i]->jacIndexes[it];
                    dep[e] = handler.createCG(Argument<Base>(*loopEnd));
                }
            }
        }

        /**
         * for each loop:
         *  ĺoop start/end
         * jac operations (1st iteration) become arguments of the loop end
         * put new operations in jac values -> an alias for the end of the loop
         * 
         * replace loop results with an actual result from an atomic loop evaluation
         */
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
    OperationNode<Base>* CLangCompileModelHelper<Base>::findSparseJacLoopResult(OperationNode<Base>* jacNode,
                                                                                std::map<OperationNode<Base>*, std::map<size_t, Argument<Base> > >& pxArgs) {
        // this is only used for the reverse mode
        if (jacNode == NULL)
            return NULL;

        const std::vector<Argument<Base> >& args = jacNode->getArguments();
        if (jacNode->getOperationType() == CGAddOp) {
            OperationNode<Base>* loop1 = findSparseJacLoopResult(args[0].getOperation(), pxArgs);
            OperationNode<Base>* loop2 = findSparseJacLoopResult(args[1].getOperation(), pxArgs);
            assert(loop1 == NULL || loop2 == NULL || loop1 == loop2);
            if (loop1 != NULL)
                return loop1;
            else
                return loop2;
        }

        OperationNode<Base>* resultNode;
        Argument<Base> argPx;
        if (jacNode->getOperationType() == CGMulOp) {
            // TODO: maybe also consider args[1]
            resultNode = args[0].getOperation();
            argPx = args[1];
        } else if (jacNode->getOperationType() == CGLoopAtomicResultOp) {
            resultNode = jacNode;
            argPx = Argument<Base>(Base(1.0));
        } else {
            return NULL;
        }

        if (resultNode != NULL && resultNode->getOperationType() == CGLoopAtomicResultOp) {
            size_t atomicJ = resultNode->getInfo()[0]; // temporary variable
            OperationNode<Base>* loopNode = resultNode->getArguments()[0].getOperation();

            typename std::map<OperationNode<Base>*, std::map<size_t, Argument<Base> > >::iterator itPxs;
            itPxs = pxArgs.find(loopNode);
            if (itPxs == pxArgs.end()) {
                pxArgs[loopNode][atomicJ] = argPx;
            } else {
                std::map<size_t, Argument<Base> >& pxs = itPxs->second;
                typename std::map<size_t, Argument<Base> > ::iterator itpx = pxs.find(atomicJ);
                if (itpx == pxs.end()) {
                    pxs[atomicJ] = argPx;
                } else {
                    OperationNode<Base>* addNode = new OperationNode<Base>(CGAddOp, itpx->second, argPx); /////////// manager.handleNode()
                    itpx->second = Argument<Base>(*addNode);
                }
            }

            return loopNode;
        }

        return NULL;
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::evalLoopTape(CodeHandler<Base>& handler,
                                                                  LoopAtomicFun<Base>& atomic,
                                                                  const OperationNode<Base>& loopEvalNode,
                                                                  const vector<OperationNode<Base>*>& indexedIndependents) {
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
                return generateForward1Graph(handler, atomic, indexedIndependents, loopEvalNode.getArguments());
            }
        } else {
            /**
             * Reverse mode
             */
            if (p == 0) {
                return generateReverse1Graph(handler, atomic, indexedIndependents, loopEvalNode.getArguments());
            }
        }
        assert(false); // TODO
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
                                                                           const std::vector<Argument<Base> >& argsAtomic) {//argsAtomic = [txAtomic]
        typedef CppAD::CG<Base> CGB;

        //size_t m = atomic.getTapeDependentCount();
        const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = atomic.getIndexedIndepIndexes();
        const std::vector<LoopPosition>& nonIndexedIndepIndexes = atomic.getNonIndexedIndepIndexes();
        const std::vector<LoopPositionTmp>& temporaryIndependents = atomic.getTemporaryIndependents();


        ADFun<CGB>* fun = atomic.getTape();

        //size_t mTape = fun->Range();
        size_t nTape = fun->Domain();
        size_t nIndexed = indexedIndepIndexes.size();
        size_t nNonIndexed = nonIndexedIndepIndexes.size();
        size_t nTmp = temporaryIndependents.size();

        // zero order
        vector<CGB> x = createLoopGraphIndependentVector(handler, atomic, indexedIndependents, argsAtomic, 1);
        assert(x.size() == nTape);

        fun->Forward(0, x);

        // forward first order
        vector<CGB> tx(nTape * 2);
        for (size_t j = 0; j < nTape; j++) {
            tx[j * 2] = x[j];
        }

        for (size_t j = 0; j < nIndexed; j++) {
            const LoopPosition& pos = indexedIndepIndexes[j][0];
            tx[pos.tape * 2 + 1] = handler.createCG(argsAtomic[pos.atomic * 2 + 1]);
        }
        for (size_t j = 0; j < nNonIndexed; j++) {
            const LoopPosition& pos = nonIndexedIndepIndexes[j];
            tx[pos.tape * 2 + 1] = handler.createCG(argsAtomic[pos.atomic * 2 + 1]);
        }
        for (size_t j = 0; j < nTmp; j++) {
            const LoopPositionTmp& pos = temporaryIndependents[j];
            tx[pos.tape * 2 + 1] = handler.createCG(argsAtomic[pos.atomic * 2 + 1]);
        }

        vector<CGB> ty = fun->Forward(1, tx);

        return ty;
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::generateReverse1Graph(CodeHandler<Base>& handler,
                                                                           LoopAtomicFun<Base>& atomic,
                                                                           const vector<OperationNode<Base>*>& indexedIndependents,
                                                                           const std::vector<Argument<Base> >& argsAtomic) { //argsAtomic = [txAtomic pyAtomic]
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
                info[0] = 1;
                info[1] = i;
                std::vector<Argument<Base> > emptyArgs;
                OperationNode<Base>* indexedPy = new OperationNode<Base>(CGLoopIndexedIndepOp, info, emptyArgs);
                py[i] = handler.createCG(indexedPy);
            }
        }

        vector<CGB> px = fun->Reverse(1, py);

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

}

#endif
