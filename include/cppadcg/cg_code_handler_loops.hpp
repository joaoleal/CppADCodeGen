#ifndef CPPAD_CG_CODE_HANDLER_LOOPS_INCLUDED
#define CPPAD_CG_CODE_HANDLER_LOOPS_INCLUDED
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

    template<class Base>
    class JacOrigElementLoopInfo {
    public:
        std::vector<size_t> jacIndexes;
        IndexPattern* jacPattern;

        inline JacOrigElementLoopInfo() :
            jacPattern(NULL) {
        }
    };

    template<class Base>
    class JacOrigElementIndepLoopInfo {
    public:
        Argument<Base> arg; // tx1 for forward, jac element for reverse
        JacOrigElementLoopInfo<Base>* origJacElement; // original model jacobian element

        inline JacOrigElementIndepLoopInfo() :
            origJacElement(NULL) {
        }
    };

    template<class Base>
    class JacTapeElementLoopInfo {
    public:
        std::map<size_t, JacOrigElementIndepLoopInfo<Base> > origIndep2Info;
    };

    template<class Base>
    vector<CG<Base> > CodeHandler<Base>::createLoopGraphIndependentVector(LoopAtomicFun<Base>& atomic,
                                                                          CodeHandler<Base>::LoopOperationGraph& graph,
                                                                          const std::vector<Argument<Base> >& args,
                                                                          size_t p) {

        typedef CppAD::CG<Base> CGB;

        const std::vector<LoopPositionTmp>& temporaryIndependents = atomic.getTemporaryIndependents();
        const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = atomic.getIndexedIndepIndexes();
        const std::vector<LoopPosition>& nonIndexedIndepIndexes = atomic.getNonIndexedIndepIndexes();

        ADFun<CGB>* fun = atomic.getTape();

        //size_t mTape = fun_->Range();
        size_t nTape = fun->Domain();
        size_t nIndexed = indexedIndepIndexes.size();
        size_t nNonIndexed = nonIndexedIndepIndexes.size();

        // indexed independents
        vector<CGB> x(nTape);
        for (size_t j = 0; j < nIndexed; j++) {
            x[j] = CGB(*this, Argument<Base>(*graph.indexedIndependents[j]));
        }
        // non indexed
        for (size_t j = 0; j < nNonIndexed; j++) {
            x[nIndexed + j] = CGB(*this, args[nonIndexedIndepIndexes[j].atomic * (p + 1)]);
        }
        // temporaries
        for (size_t j = 0; j < temporaryIndependents.size(); j++) {
            x[nIndexed + nNonIndexed + j] = CGB(*this, args[temporaryIndependents[j].atomic * (p + 1)]);
        }
        return x;
    }

    template<class Base>
    void CodeHandler<Base>::registerLoop(LoopAtomicFun<Base>& loop) {
        typename std::map<size_t, LoopData*>::const_iterator it = _loops.find(loop.getLoopId());
        if (it != _loops.end()) {
            delete it->second;
        }
        _loops[loop.getLoopId()] = new LoopData(loop);
    }

    template<class Base>
    LoopAtomicFun<Base>* CodeHandler<Base>::getLoop(size_t loopId) const {
        typename std::map<size_t, LoopData*>::const_iterator it = _loops.find(loopId);
        if (it != _loops.end()) {
            return it->second->atomic;
        }

        return NULL;
    }

    template<class Base>
    void CodeHandler<Base>::prepareLoops(std::vector<CG<Base> >& dependent) {

        std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > > jacIndexPatterns;
        std::map<OperationNode<Base>*, OperationNode<Base>*> loopDependents;
        for (size_t i = 0; i < dependent.size(); i++) {
            if (dependent[i].getOperationNode() != NULL) {
                OperationNode<Base>* result = dependent[i].getOperationNode();
                if (result->getOperationType() == CGLoopResultOp) {
                    OperationNode<Base>* loop = result->getArguments()[0].getOperation();
                    loopDependents[result] = loop;
                }
            }
        }

        prepareLoops(loopDependents, jacIndexPatterns);
    }

    template<class Base>
    void CodeHandler<Base>::prepareLoops(const std::map<OperationNode<Base>*, OperationNode<Base>*>& loopDependents,
                                         const std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > >& jacIndexPatterns) {
        typename std::map<OperationNode<Base>*, OperationNode<Base>*>::const_iterator it;
        for (it = loopDependents.begin(); it != loopDependents.end(); ++it) {
            OperationNode<Base>* result = it->first;
            OperationNode<Base>* loop = it->second;

            insertLoopOperations(*result, *loop, jacIndexPatterns);
        }
    }

    template<class Base>
    void CodeHandler<Base>::insertLoopOperations(OperationNode<Base>& loopResult,
                                                 OperationNode<Base>& loop,
                                                 const std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > >& jacIndexPatterns) {
        assert(loop.getOperationType() == CGLoopForwardOp || loop.getOperationType() == CGLoopReverseOp);
        assert(loop.getInfo().size() == 4);

        size_t loopId = loop.getInfo()[0];
        size_t p = loop.getInfo()[3];

        loopResult.operation_ = CGLoopResultOp;
        loopResult.arguments_.clear();
        loopResult.info_.clear();
        LoopData* loopData = _loops.at(loopId);

        if (loop.getOperationType() == CGLoopForwardOp) {
            /**
             * Forward mode
             */
            if (p == 0) {
                if (loopData->graphForward0 == NULL) {
                    loopData->graphForward0 = generateLoopForward0Graph(*loopData->atomic, p, loop.getArguments());
                    assert(loopData->graphForward0 != NULL);
                }

                loopResult.arguments_.push_back(Argument<Base>(*loopData->graphForward0->loopEnd));

            } else if (p == 1) {
                generateForward1Graph(*loopData, jacIndexPatterns.at(loopData->atomic), loop.getArguments(), false); //////////////
                assert(loopData->graphForward1 != NULL);

                loopResult.arguments_.push_back(Argument<Base>(*loopData->graphForward1->loopEnd));
            } else {
                assert(false); // TODO
            }
        } else {
            /**
             * Reverse mode
             */
            if (p == 0) {
                generateReverse1Graph(*loopData, jacIndexPatterns.at(loopData->atomic), loop.getArguments(), false); //////////////
                assert(loopData->graphReverse1 != NULL);

                loopResult.arguments_.push_back(Argument<Base>(*loopData->graphReverse1->loopEnd));
            } else {
                assert(false); // TODO
            }
        }

    }

    template<class Base>
    typename CodeHandler<Base>::LoopOperationGraph* CodeHandler<Base>::generateLoopForward0Graph(LoopAtomicFun<Base>& atomic,
                                                                                                 size_t p,
                                                                                                 const std::vector<Argument<Base> >& args) {
        typedef CppAD::CG<Base> CGB;

        size_t iterationCount = atomic.getIterationCount();
        const std::vector<LoopPositionTmp>& temporaryIndependents = atomic.getTemporaryIndependents();
        const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = atomic.getIndexedIndepIndexes();
        //const std::vector<LoopPosition>& nonIndexedIndepIndexes = atomic.getNonIndexedIndepIndexes();
        ADFun<CGB>*fun = atomic.getTape();

        //size_t mTape = fun_->Range();
        //size_t nTape = fun->Domain();
        std::vector<size_t> startEndInfo(2);
        startEndInfo[0] = atomic.getLoopId();
        startEndInfo[1] = iterationCount;

        std::vector<size_t> info(1);

        LoopOperationGraph* graphForward = new LoopOperationGraph();
        std::vector<Argument<Base> > startArgs(temporaryIndependents.size());
        for (size_t j = 0; j < temporaryIndependents.size(); j++) {
            const LoopPositionTmp& pos = temporaryIndependents[j];
            startArgs[j] = args[pos.atomic];
        }

        graphForward->loopStart = new OperationNode<Base>(CGLoopStartOp, startEndInfo, startArgs);
        manageOperationNode(graphForward->loopStart);

        // indexed independents
        size_t nIndexed = indexedIndepIndexes.size();
        graphForward->indexedIndependents.resize(nIndexed);

        std::vector<Argument<Base> > indexedArgs(iterationCount + 1);
        indexedArgs[0] = Argument<Base>(*graphForward->loopStart);
        for (size_t j = 0; j < nIndexed; j++) {
            for (size_t it = 0; it < iterationCount; it++) {
                //assert(args[j].getOperation() != NULL);
                indexedArgs[it + 1] = args[indexedIndepIndexes[j][it].atomic];
            }
            info[0] = j;
            graphForward->indexedIndependents[j] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, indexedArgs);
            manageOperationNode(graphForward->indexedIndependents[j]);
        }

        vector<CGB> tx = createLoopGraphIndependentVector(atomic, *graphForward, args, 0);
        vector<CGB> ty = fun->Forward(p, tx);

        size_t ty_size = ty.size();
        graphForward->indexedResults.resize(ty_size);
        std::vector<Argument<Base> > endArgs(ty_size);
        indexedArgs.resize(1);
        for (size_t i = 0; i < ty_size; i++) {
            indexedArgs[0] = asArgument(ty[i]);
            info[0] = _loopDependentIndexPatterns.size(); // dependent index pattern location
            _loopDependentIndexPatterns.push_back(atomic.getDependentIndexPatterns()[i]);
            OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, info, indexedArgs);
            graphForward->indexedResults[i] = yIndexed;
            manageOperationNode(yIndexed);
            endArgs[i] = Argument<Base>(*yIndexed);
        }

        graphForward->loopEnd = new OperationNode<Base>(CGLoopEndOp, startEndInfo, endArgs);
        manageOperationNode(graphForward->loopEnd);

        return graphForward;
    }

    template<class Base>
    void CodeHandler<Base>::generateForward1Graph(LoopData& loopData,
                                                  const std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > >& jac,
                                                  const std::vector<Argument<Base> >& args,
                                                  bool variableTx1) {
        LoopAtomicFun<Base>* atomic = loopData.atomic;
        size_t iterationCount = atomic->getIterationCount();
        //size_t mFull = atomic->getLoopDependentCount();
        //size_t nFull = atomic->getLoopIndependentCount();
        //size_t m = atomic->getTapeDependentCount();
        //size_t n = atomic->getTapeIndependentCount();
        const std::vector<std::vector<LoopPosition> >& indexIndep = atomic->getIndexedIndepIndexes();
        const std::vector<LoopPosition>& nonindexIndep = atomic->getNonIndexedIndepIndexes();
        const std::vector<LoopPositionTmp>& tempIndep = atomic->getTemporaryIndependents();

        // transpose jac
        // loop -> tape independent -> equation pattern -> orig independent(temporaries only) -> iteration = position
        std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > jacT;
        typename std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > >::const_iterator itr;
        for (itr = jac.begin(); itr != jac.end(); ++itr) {
            size_t r = itr->first;
            typename std::map<size_t, JacTapeElementLoopInfo<Base> >::const_iterator itc;
            for (itc = itr->second.begin(); itc != itr->second.end(); ++itc) {
                size_t c = itc->first;
                jacT[c][r] = itc->second;
            }
        }

        /**
         * args = [ tx ]
         */
        assert(args.size() == 2 * atomic->getLoopIndependentCount());

        std::map<const IndexPattern*, CG<Base> > jacVals;

        // indexed
        for (size_t jTape = 0; jTape < indexIndep.size(); jTape++) {
            for (size_t it = 0; it < iterationCount; it++) {
                const LoopPosition& pos = indexIndep[jTape][it];
                const Argument<Base>& argTx1 = args[pos.atomic * 2 + 1];
                if (argTx1.getParameter() == NULL || !IdenticalZero(*argTx1.getParameter())) {
                    if (loopData.graphForward1Indeps.find(jTape) == loopData.graphForward1Indeps.end()) {
                        loopData.graphForward1 = generateForward1Graph(*atomic, loopData.graphForward1, jTape, jacT.at(jTape), jacVals, args, variableTx1);
                        loopData.graphForward1Indeps.insert(jTape);
                    }
                    break;
                }
            }
        }

        // non-indexed
        for (size_t j = 0; j < nonindexIndep.size(); j++) {
            const LoopPosition& pos = nonindexIndep[j];
            size_t jTape = indexIndep.size() + j;
            const Argument<Base>& argTx1 = args[pos.atomic * 2 + 1];
            if (argTx1.getParameter() == NULL || !IdenticalZero(*argTx1.getParameter())) {
                if (loopData.graphForward1Indeps.find(jTape) == loopData.graphForward1Indeps.end()) {
                    loopData.graphForward1 = generateForward1Graph(*atomic, loopData.graphForward1, jTape, jacT.at(jTape), jacVals, args, variableTx1);
                    loopData.graphForward1Indeps.insert(jTape);
                }
                break;
            }
        }

        // temporary
        for (size_t j = 0; j < tempIndep.size(); j++) {
            const LoopPositionTmp& pos = tempIndep[j];
            size_t jTape = indexIndep.size() + nonindexIndep.size() + j;
            const Argument<Base>& argTx1 = args[pos.atomic * 2 + 1];
            if (argTx1.getParameter() == NULL || !IdenticalZero(*argTx1.getParameter())) {
                if (loopData.graphForward1Indeps.find(jTape) == loopData.graphForward1Indeps.end()) {
                    loopData.graphForward1 = generateForward1Graph(*atomic, loopData.graphForward1, jTape, jacT.at(jTape), jacVals, args, variableTx1);
                    loopData.graphForward1Indeps.insert(jTape);
                }
                break;
            }
        }

        assert(loopData.graphForward1 != NULL);

        std::vector<Argument<Base> > tyIndexedArgs(1);
        std::vector<size_t> info(1);

        typename std::map<const IndexPattern*, CG<Base> >::const_iterator itv;
        for (itv = jacVals.begin(); itv != jacVals.end(); ++itv) {
            const IndexPattern* ip = itv->first;

            tyIndexedArgs[0] = asArgument(itv->second);

            info[0] = _loopDependentIndexPatterns.size(); // dependent index pattern location
            _loopDependentIndexPatterns.push_back(ip);

            OperationNode<Base>* tyIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, info, tyIndexedArgs);
            manageOperationNode(tyIndexed);

            loopData.graphForward1->indexedResults.push_back(tyIndexed); // is this really needed?????????????????
            loopData.graphForward1->loopEnd->arguments_.push_back(Argument<Base>(*tyIndexed));
        }

    }

    template<class Base>
    typename CodeHandler<Base>::LoopOperationGraph* CodeHandler<Base>::generateForward1Graph(LoopAtomicFun<Base>& atomic,
                                                                                             LoopOperationGraph* graphForward1,
                                                                                             size_t jTape,
                                                                                             const std::map<size_t, JacTapeElementLoopInfo<Base> >& jacCol,
                                                                                             std::map<const IndexPattern*, CG<Base> >& jacColVals,
                                                                                             const std::vector<Argument<Base> >& args,
                                                                                             bool variableTx1) {
        typedef CppAD::CG<Base> CGB;

        size_t iterationCount = atomic.getIterationCount();
        //size_t m = atomic.getTapeDependentCount();
        //const std::vector<LoopPositionTmp>& temporaryIndependents = atomic.getTemporaryIndependents();
        const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = atomic.getIndexedIndepIndexes();
        //const std::vector<LoopPosition>& nonIndexedIndepIndexes = atomic.getNonIndexedIndepIndexes();

        ADFun<CGB>*fun = atomic.getTape();

        //size_t mTape = fun->Range();
        size_t nTape = fun->Domain();
        size_t nIndexed = indexedIndepIndexes.size();
        //size_t nNonIndexed = nonIndexedIndepIndexes.size();

        if (graphForward1 == NULL) {
            graphForward1 = new LoopOperationGraph();

            std::vector<size_t> startEndInfo(2);
            startEndInfo[0] = atomic.getLoopId();
            startEndInfo[1] = iterationCount;

            std::vector<Argument<Base> > startArgs; // will be updated later (temporary variables)

            graphForward1->loopStart = new OperationNode<Base>(CGLoopStartOp, startEndInfo, startArgs);
            manageOperationNode(graphForward1->loopStart);

            std::vector<size_t> info(1);

            // indexed independents
            graphForward1->indexedIndependents.resize(nIndexed);
            std::vector<Argument<Base> > xIndexedArgs(iterationCount + 1);
            xIndexedArgs[0] = Argument<Base>(*graphForward1->loopStart);
            for (size_t j = 0; j < nIndexed; j++) {
                for (size_t it = 0; it < iterationCount; it++) {
                    xIndexedArgs[it + 1] = args[indexedIndepIndexes[j][it].atomic];
                }
                info[0] = j;
                graphForward1->indexedIndependents[j] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, xIndexedArgs);
                manageOperationNode(graphForward1->indexedIndependents[j]);
            }

            graphForward1->indexedTx1.resize(nTape);
            info.resize(2);
            info[0] = 1; // how to know its tx1 and not x
            std::vector<Argument<Base> > emptyArgs;
            for (size_t j = 0; j < nTape; j++) {
                info[1] = j;
                graphForward1->indexedTx1[j] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, emptyArgs);
                manageOperationNode(graphForward1->indexedTx1[j]);
            }

            std::vector<Argument<Base> > endArgs; // will be updated later with the results
            graphForward1->loopEnd = new OperationNode<Base>(CGLoopEndOp, startEndInfo, endArgs);
            manageOperationNode(graphForward1->loopEnd);
        }

        // zero order
        vector<CGB> x = createLoopGraphIndependentVector(atomic, *graphForward1, args, 1);
        assert(x.size() == nTape);

        fun->Forward(0, x);

        // forward first order
        vector<CGB> tx(nTape * 2);
        for (size_t j = 0; j < nTape; j++) {
            tx[j * 2] = x[j];
            tx[j * 2 + 1] = Base(0);
        }

        if (variableTx1)
            tx[jTape * 2 + 1] = CGB(*this, Argument<Base>(*graphForward1->indexedTx1[jTape]));
        else
            tx[jTape * 2 + 1] = Base(1);
        vector<CGB> ty = fun->Forward(1, tx);

        graphForward1->indexedResults.reserve(graphForward1->indexedResults.size() + jacCol.size());
        graphForward1->loopEnd->arguments_.reserve(graphForward1->loopEnd->arguments_.size() + jacCol.size());

        typename std::map<size_t, JacTapeElementLoopInfo<Base> >::const_iterator itJac;
        for (itJac = jacCol.begin(); itJac != jacCol.end(); ++itJac) {
            size_t tapeJ = itJac->first;
            const JacTapeElementLoopInfo<Base>& jacTapeInfo = itJac->second;

            typename std::map<size_t, JacOrigElementIndepLoopInfo<Base> >::const_iterator itO;
            for (itO = jacTapeInfo.origIndep2Info.begin(); itO != jacTapeInfo.origIndep2Info.end(); ++itO) {
                const JacOrigElementIndepLoopInfo<Base>& origInfo = itO->second;

                if (origInfo.arg.getOperation() != NULL && origInfo.arg.getOperation()->getOperationType() != CGInvOp) {
                    graphForward1->loopStart->arguments_.push_back(origInfo.arg);
                }

                const IndexPattern* ip = origInfo.origJacElement->jacPattern;

                CGB tx = CGB(*this, origInfo.arg);
                typename std::map<const IndexPattern*, CGB>::iterator itJacEl = jacColVals.find(ip);
                if (itJacEl == jacColVals.end()) {
                    jacColVals[ip] = ty[tapeJ * 2 + 1] * tx;
                } else {
                    itJacEl->second += ty[tapeJ * 2 + 1] * tx;
                }
            }
        }

        /**
         * TODO:
         * 
         * must know if the temporary variables were needed!
         * If so then they must be made a dependency of the loop start
         */

        /**
         * Non-indexed operations should go outside the loop!
         */

        return graphForward1;
    }

    template<class Base>
    void CodeHandler<Base>::generateReverse1Graph(LoopData& loopData,
                                                  const std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > >& jac,
                                                  const std::vector<Argument<Base> >& args,
                                                  bool variablePy) {
        LoopAtomicFun<Base>* atomic = loopData.atomic;
        size_t iterationCount = atomic->getIterationCount();
        //size_t mFull = atomic->getLoopDependentCount();
        size_t nFull = atomic->getLoopIndependentCount();
        size_t m = atomic->getTapeDependentCount();
        const std::vector<std::vector<LoopPosition> >& dependentIndexes = atomic->getDependentIndexes();

        /**
         * args = [ x , py ]
         */
        assert(args.size() == nFull + atomic->getLoopDependentCount());

        for (size_t i = 0; i < m; i++) {
            for (size_t it = 0; it < iterationCount; it++) {
                const Argument<Base>& argPy = args[nFull + dependentIndexes[i][it].atomic];
                if (argPy.getParameter() == NULL || !IdenticalZero(*argPy.getParameter())) {
                    if (loopData.graphReverse1Deps.find(i) == loopData.graphReverse1Deps.end()) {
                        loopData.graphReverse1 = generateReverse1Graph(*atomic, loopData.graphReverse1, i, jac.at(i), args, variablePy);
                        loopData.graphReverse1Deps.insert(i);
                    }
                    break;
                }
            }
        }
    }

    template<class Base>
    typename CodeHandler<Base>::LoopOperationGraph* CodeHandler<Base>::generateReverse1Graph(LoopAtomicFun<Base>& atomic,
                                                                                             LoopOperationGraph* graphReverse1,
                                                                                             size_t dep,
                                                                                             const std::map<size_t, JacTapeElementLoopInfo<Base> >& jacRow,
                                                                                             const std::vector<Argument<Base> >& args,
                                                                                             bool variablePy) {
        typedef CppAD::CG<Base> CGB;

        std::map<const IndexPattern*, CG<Base> > jacRowVals;

        size_t iterationCount = atomic.getIterationCount();
        size_t m = atomic.getTapeDependentCount();
        //const std::vector<LoopPositionTmp>& temporaryIndependents = atomic.getTemporaryIndependents();
        const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = atomic.getIndexedIndepIndexes();
        //const std::vector<LoopPosition>& nonIndexedIndepIndexes = atomic.getNonIndexedIndepIndexes();

        ADFun<CGB>*fun = atomic.getTape();

        //size_t mTape = fun_->Range();
        //size_t nTape = fun->Domain();
        size_t nIndexed = indexedIndepIndexes.size();
        //size_t nNonIndexed = nonIndexedIndepIndexes.size();

        if (graphReverse1 == NULL) {
            graphReverse1 = new LoopOperationGraph();

            std::vector<size_t> startEndInfo(2);
            startEndInfo[0] = atomic.getLoopId();
            startEndInfo[1] = iterationCount;

            std::vector<Argument<Base> > startArgs; // will be updated later (temporary variables)

            graphReverse1->loopStart = new OperationNode<Base>(CGLoopStartOp, startEndInfo, startArgs);
            manageOperationNode(graphReverse1->loopStart);

            std::vector<size_t> info(1);

            // indexed independents
            graphReverse1->indexedIndependents.resize(nIndexed);
            std::vector<Argument<Base> > xIndexedArgs(iterationCount + 1);
            xIndexedArgs[0] = Argument<Base>(*graphReverse1->loopStart);
            for (size_t j = 0; j < nIndexed; j++) {
                for (size_t it = 0; it < iterationCount; it++) {
                    xIndexedArgs[it + 1] = args[indexedIndepIndexes[j][it].atomic];
                }
                info[0] = j;
                graphReverse1->indexedIndependents[j] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, xIndexedArgs);
                manageOperationNode(graphReverse1->indexedIndependents[j]);
            }

            graphReverse1->indexedPy.resize(m);
            info.resize(2);
            info[0] = 1; // how to know its py and not x
            std::vector<Argument<Base> > emptyArgs;
            for (size_t i = 0; i < m; i++) {
                info[1] = i;
                graphReverse1->indexedPy[i] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, emptyArgs);
                manageOperationNode(graphReverse1->indexedPy[i]);
            }

            std::vector<Argument<Base> > endArgs; // will be updated later with the results
            graphReverse1->loopEnd = new OperationNode<Base>(CGLoopEndOp, startEndInfo, endArgs);
            manageOperationNode(graphReverse1->loopEnd);
        }

        // zero order
        vector<CGB> x = createLoopGraphIndependentVector(atomic, *graphReverse1, args, 0);
        assert(x.size() == fun->Domain());

        fun->Forward(0, x);

        // reverse first order
        vector<CGB> py(m);
        if (variablePy)
            py[dep] = CGB(*this, Argument<Base>(*graphReverse1->indexedPy[dep]));
        else
            py[dep] = Base(1);
        vector<CGB> px = fun->Reverse(1, py);

        graphReverse1->indexedResults.reserve(graphReverse1->indexedResults.size() + jacRow.size());
        graphReverse1->loopEnd->arguments_.reserve(graphReverse1->loopEnd->arguments_.size() + jacRow.size());

        typename std::map<size_t, JacTapeElementLoopInfo<Base> >::const_iterator itJac;
        for (itJac = jacRow.begin(); itJac != jacRow.end(); ++itJac) {
            size_t j = itJac->first;
            const JacTapeElementLoopInfo<Base>& jacTapeInfo = itJac->second;

            typename std::map<size_t, JacOrigElementIndepLoopInfo<Base> >::const_iterator itO;
            for (itO = jacTapeInfo.origIndep2Info.begin(); itO != jacTapeInfo.origIndep2Info.end(); ++itO) {
                const JacOrigElementIndepLoopInfo<Base>& origInfo = itO->second;

                if (origInfo.arg.getOperation() != NULL && origInfo.arg.getOperation()->getOperationType() != CGInvOp) {
                    graphReverse1->loopStart->arguments_.push_back(origInfo.arg);
                }

                const IndexPattern* ip = origInfo.origJacElement->jacPattern;

                CGB py = CGB(*this, origInfo.arg);
                typename std::map<const IndexPattern*, CGB>::iterator itJacEl = jacRowVals.find(ip);
                if (itJacEl == jacRowVals.end()) {
                    jacRowVals[ip] = px[j] * py;
                } else {
                    itJacEl->second += px[j] * py;
                }
            }
        }

        std::vector<Argument<Base> > pxIndexedArgs(1);
        std::vector<size_t> info(1);

        typename std::map<const IndexPattern*, CG<Base> >::const_iterator itv;
        for (itv = jacRowVals.begin(); itv != jacRowVals.end(); ++itv) {
            const IndexPattern* ip = itv->first;

            pxIndexedArgs[0] = asArgument(itv->second);

            info[0] = _loopDependentIndexPatterns.size(); // dependent index pattern location
            _loopDependentIndexPatterns.push_back(ip);

            OperationNode<Base>* pxIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, info, pxIndexedArgs);
            manageOperationNode(pxIndexed);

            graphReverse1->indexedResults.push_back(pxIndexed); // is this really needed?????????????????
            graphReverse1->loopEnd->arguments_.push_back(Argument<Base>(*pxIndexed));
        }

        /**
         * TODO:
         * 
         * must know if the temporary variables were needed!
         * If so then they must be made a dependency of the loop start
         */

        /**
         * Non-indexed operations should go outside the loop!
         */

        return graphReverse1;
    }

}

#endif