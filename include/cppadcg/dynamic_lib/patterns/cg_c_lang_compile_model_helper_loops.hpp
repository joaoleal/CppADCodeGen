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

    /***************************************************************************
     *  Methods related with loop insertion into the operation graph
     **************************************************************************/

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareLoops(CodeHandler<Base>& handler,
                                                     std::vector<CGBase>& dep,
                                                     std::map<LoopAtomicFun<Base>*, std::map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> > >& evaluations1it,
                                                     std::map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo<Base>* > >& dependentIndexes,
                                                     size_t assignOrAdd) {

        std::map<LoopAtomicFun<Base>*, vector<OperationNode<Base>* > > loopIndexedIndependents;
        std::map<LoopAtomicFun<Base>*, vector<OperationNode<Base>* > > loopIndexedIndependents2;

        /**
         * insert the loops into the operation graph
         */
        typename std::map<LoopAtomicFun<Base>*, std::map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> > >::iterator itL;
        for (itL = evaluations1it.begin(); itL != evaluations1it.end(); ++itL) {
            LoopAtomicFun<Base>* loopFunc = itL->first;
            std::map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> >& loopAtomEvaluations1it = itL->second;

            /**
             * make the loop start
             */
            OperationNode<Base>* loopStart = new LoopStartOperationNode<Base>(*loopFunc);
            handler.manageOperationNodeMemory(loopStart);

            IndexOperationNode<Base>* iterationIndexOp = new IndexOperationNode<Base>(LoopAtomicFun<Base>::ITERATION_INDEX, *loopStart);
            handler.manageOperationNodeMemory(iterationIndexOp);
            std::set<IndexOperationNode<Base>*> indexesOps;
            indexesOps.insert(iterationIndexOp);

            // indexed independents
            vector<OperationNode<Base>* >& indexedIndependents = loopIndexedIndependents[loopFunc]; // zero order
            vector<OperationNode<Base>* >& indexedIndependents2 = loopIndexedIndependents2[loopFunc]; // first order

            /**
             * generate the expressions inside the loop
             */
            replaceAtomicLoopWithExpression(handler, *loopFunc,
                                            *iterationIndexOp, loopAtomEvaluations1it,
                                            indexedIndependents, indexedIndependents2);
            /**
             * make the loop end
             */
            const vector<IndexedDependentLoopInfo<Base>* >& dependents = dependentIndexes.at(loopFunc);
            size_t dep_size = dependents.size();
            vector<std::pair<CGBase, IndexPattern*> > indexedLoopResults(dep_size);
            for (size_t i = 0; i < dep_size; i++) {
                indexedLoopResults[i] = std::make_pair(dependents[i]->origVals[0],
                                                       dependents[i]->pattern);
            }
            OperationNode<Base>* loopEnd = createLoopEnd(handler, indexedLoopResults, indexesOps, *loopFunc, assignOrAdd);

            /**
             * change the dependents (must depend directly on the loop)
             */
            for (size_t i = 0; i < dep_size; i++) {
                for (size_t it = 0; it < dependents[i]->indexes.size(); it++) {
                    size_t e = dependents[i]->indexes[it];
                    dep[e] = handler.createCG(Argument<Base>(*loopEnd));
                }
            }

            /**
             * move no-nindexed expressions outside loop
             */
            moveNonIndexedOutsideLoop(*loopStart, *loopEnd, LoopAtomicFun<Base>::ITERATION_INDEX);
        }

    }

    template<class Base>
    void CLangCompileModelHelper<Base>::replaceAtomicLoopWithExpression(CodeHandler<Base>& handler,
                                                                        LoopAtomicFun<Base>& loopFunc,
                                                                        IndexOperationNode<Base>& iterationIndexOp,
                                                                        std::map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> >& evaluations1it,
                                                                        vector<OperationNode<Base>* >& indexedIndependents,
                                                                        vector<OperationNode<Base>* >& indexedIndependents2) {
        /**
         * insert the loops into the operation graph
         */

        // indexed independents
        if (indexedIndependents.size() == 0) {
            size_t nIndexed = loopFunc.getIndependentIndexPatterns().size();
            indexedIndependents.resize(nIndexed); // zero order
            std::vector<Argument<Base> > xIndexedArgs(1);
            xIndexedArgs[0] = Argument<Base>(iterationIndexOp);
            std::vector<size_t> info(2);
            info[0] = 0; // tx
            for (size_t j = 0; j < nIndexed; j++) {
                info[1] = handler.addLoopIndependentIndexPattern(*loopFunc.getIndependentIndexPatterns()[j], j);
                indexedIndependents[j] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, xIndexedArgs);
                handler.manageOperationNodeMemory(indexedIndependents[j]);
            }
        }

        typename std::map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> >::iterator itE;
        for (itE = evaluations1it.begin(); itE != evaluations1it.end(); ++itE) {
            LoopEvaluationOperationNode<Base>* loopEvalNode = itE->first;
            vector<OperationNode<Base>*>& resultNodes1it = itE->second;

            // evaluate tape 1st iteration (generate the operation graph)
            if (loopEvalNode != NULL) {
                vector<CG<Base> > newTapeResults1it = evalLoopTape(handler, *loopEvalNode,
                                                                   indexedIndependents, indexedIndependents2,
                                                                   iterationIndexOp);

                assert(resultNodes1it.size() == newTapeResults1it.size());
                for (size_t j = 0; j < resultNodes1it.size(); j++) {
                    if (resultNodes1it[j] != NULL) {
                        assert(resultNodes1it[j]->getOperationType() == CGLoopAtomicResultOp);
                        resultNodes1it[j]->getArguments()[0] = asArgument(newTapeResults1it[j]);
                    }
                }
            }
        }

    }

    template<class Base>
    OperationNode<Base>* CLangCompileModelHelper<Base>::createLoopEnd(CodeHandler<Base>& handler,
                                                                      const vector<std::pair<CG<Base>, IndexPattern*> >& indexedLoopResults,
                                                                      const std::set<IndexOperationNode<Base>*>& indexesOps,
                                                                      const LoopNodeInfo<Base>& loopInfo,
                                                                      size_t assignOrAdd) {
        std::vector<Argument<Base> > endArgs;
        std::vector<Argument<Base> > indexedArgs(1 + indexesOps.size());
        std::vector<size_t> info(2);

        size_t dep_size = indexedLoopResults.size();
        for (size_t i = 0; i < dep_size; i++) {
            const std::pair<CG<Base>, IndexPattern*>& depInfo = indexedLoopResults[i];
            indexedArgs.resize(1);

            indexedArgs[0] = asArgument(depInfo.first); // indexed expression
            typename std::set<IndexOperationNode<Base>*>::const_iterator itIndexOp;
            for (itIndexOp = indexesOps.begin(); itIndexOp != indexesOps.end(); ++itIndexOp) {
                indexedArgs.push_back(Argument<Base>(**itIndexOp)); // dependency on the index
            }

            info[0] = handler.addLoopDependentIndexPattern(*depInfo.second); // dependent index pattern location
            info[1] = assignOrAdd;

            OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, info, indexedArgs);
            handler.manageOperationNodeMemory(yIndexed);
            endArgs.push_back(Argument<Base>(*yIndexed));
        }

        OperationNode<Base>* loopEnd = new LoopEndOperationNode<Base>(loopInfo, endArgs);
        handler.manageOperationNodeMemory(loopEnd);

        return loopEnd;
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::findLoopEvaluations(CodeHandler<Base>& handler,
                                                            OperationNode<Base>* node,
                                                            std::map<LoopEvaluationOperationNode<Base>*, std::map<size_t, OperationNode<Base>*> >& evals) {
        if (node == NULL) {
            return;
        }

        CGOpCode op = node->getOperationType();
        if (op == CGLoopAtomicResultOp) {
            size_t pos = node->getInfo()[0];
            OperationNode<Base>* loopNode = node->getArguments()[0].getOperation();
            assert(loopNode->getOperationType() == CGLoopForwardOp ||
                   loopNode->getOperationType() == CGLoopReverseOp);
            LoopEvaluationOperationNode<Base>* loopEvalNode = static_cast<LoopEvaluationOperationNode<Base>*> (loopNode);
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
                                                            std::map<LoopEvaluationOperationNode<Base>*, std::map<size_t, OperationNode<Base>*> >& evals) {
        if (node == NULL) {
            return;
        }

        CGOpCode op = node->getOperationType();
        if (op == CGLoopAtomicResultOp) {
            size_t pos = node->getInfo()[0];

            OperationNode<Base>* loopNode = node->getArguments()[0].getOperation();
            assert(loopNode->getOperationType() == CGLoopForwardOp ||
                   loopNode->getOperationType() == CGLoopReverseOp);
            LoopEvaluationOperationNode<Base>* loopEvalNode = static_cast<LoopEvaluationOperationNode<Base>*> (loopNode);

            if (loopEvalNode->getLoopAtomicFun().getLoopId() == loop->getLoopId()) {
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
                                                                  LoopEvaluationOperationNode<Base>& loopEvalNode,
                                                                  const vector<OperationNode<Base>*>& indexedIndependents,
                                                                  vector<OperationNode<Base>* >& indexedIndependents2,
                                                                  IndexOperationNode<Base>& iterationIndexOp) {
        assert(loopEvalNode.getOperationType() == CGLoopForwardOp || loopEvalNode.getOperationType() == CGLoopReverseOp);

        size_t p = loopEvalNode.getOrder();
        LoopAtomicFun<Base>& atomic = loopEvalNode.getLoopAtomicFun();

        if (loopEvalNode.getOperationType() == CGLoopForwardOp) {
            /**
             * Forward mode
             */
            if (p == 0) {
                return generateLoopForward0Graph(handler, atomic, indexedIndependents, loopEvalNode.getArguments());
            } else if (p == 1) {
                return generateForward1Graph(handler, atomic,
                                             indexedIndependents, indexedIndependents2,
                                             loopEvalNode.getArguments(), iterationIndexOp);
            }
        } else {
            /**
             * Reverse mode
             */
            if (p == 0) {
                return generateReverse1Graph(handler, atomic,
                                             indexedIndependents,
                                             loopEvalNode.getArguments(), iterationIndexOp);
            } else if (p == 1) {
                return generateReverse2Graph(handler, atomic,
                                             indexedIndependents, indexedIndependents2,
                                             loopEvalNode.getArguments(), iterationIndexOp);
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
                                                                           IndexOperationNode<Base>& iterationIndexOp) {
        typedef CppAD::CG<Base> CGB;

        ADFun<CGB>* fun = atomic.getTape();

        // zero order
        vector<CGB> x = createLoopGraphIndependentVector(handler, atomic, indexedIndependents, argsAtomic, 1);
        assert(x.size() == fun->Domain());

        fun->Forward(0, x);

        // forward first order
        vector<CGB> tx = createLoopGraphIndependentVectorTx2(handler, atomic, x, indexedIndependents2, argsAtomic, 1, iterationIndexOp);
        assert(tx.size() == fun->Domain() * 2);

        vector<CGB> ty = fun->Forward(1, tx);

        return ty;
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::generateReverse1Graph(CodeHandler<Base>& handler,
                                                                           LoopAtomicFun<Base>& atomic,
                                                                           const vector<OperationNode<Base>*>& indexedIndependents,
                                                                           const std::vector<Argument<Base> >& argsAtomic, //argsAtomic = [txAtomic pyAtomic]
                                                                           IndexOperationNode<Base>& iterationIndexOp) {
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
                info[0] = 2; // py
                info[1] = handler.addLoopIndependentIndexPattern(*atomic.getDependentIndexPatterns()[i], indexedIndependents.size() + i);
                std::vector<Argument<Base> > args(1);
                args[0] = iterationIndexOp;
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
                                                                           IndexOperationNode<Base>& iterationIndexOp) {
        typedef CppAD::CG<Base> CGB;

        size_t m = atomic.getTapeDependentCount();
        size_t nFull = atomic.getLoopIndependentCount();

        ADFun<CGB>* fun = atomic.getTape();

        // zero order
        vector<CGB> x = createLoopGraphIndependentVector(handler, atomic, indexedIndependents, argsAtomic, 1);
        assert(x.size() == fun->Domain());

        fun->Forward(0, x);

        // forward first order
        vector<CGB> tx = createLoopGraphIndependentVectorTx2(handler, atomic, x, indexedIndependents2, argsAtomic, 1, iterationIndexOp);
        assert(tx.size() == fun->Domain() * 2);

        fun->Forward(1, tx);

        // reverse second order
        vector<CGB> py(2 * m);
        for (size_t i = 0; i < 2 * m; i++) {
            if (argsAtomic[nFull * 2 + i].getOperation() == NULL) {
                py[i] = handler.createCG(argsAtomic[nFull * 2 + i]);
            } else {
                size_t i2 = i / 2; // TODO: improve this
                std::vector<size_t> info(2);
                info[0] = 2; // py
                info[1] = handler.addLoopIndependentIndexPattern(*atomic.getDependentIndexPatterns()[i2], indexedIndependents.size() + i2);
                std::vector<Argument<Base> > args(1);
                args[0] = iterationIndexOp;
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
                                                                                         IndexOperationNode<Base>& iterationIndexOp) {
        typedef CppAD::CG<Base> CGB;

        const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = atomic.getIndexedIndepIndexes();
        const std::vector<LoopPosition>& nonIndexedIndepIndexes = atomic.getNonIndexedIndepIndexes();
        const std::vector<LoopPositionTmp>& temporaryIndependents = atomic.getTemporaryIndependents();

        if (indexedIndependents2.size() == 0) {
            size_t nIndexed = indexedIndepIndexes.size();
            indexedIndependents2.resize(nIndexed); // first order
            std::vector<Argument<Base> > xIndexedArgs(1);
            xIndexedArgs[0] = Argument<Base>(iterationIndexOp);
            std::vector<size_t> info(2);
            info[0] = 1; // tx1
            for (size_t j = 0; j < nIndexed; j++) {
                info[1] = handler.addLoopIndependentIndexPattern(*atomic.getIndependentIndexPatterns()[j], j);
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
                                                                  OperationNode<Base>& loopEnd,
                                                                  const Index& loopIndex) {
        //EquationPattern<Base>::uncolor(dependents[dep].getOperationNode());
        std::set<OperationNode<Base>*> nonIndexed;

        const std::vector<Argument<Base> >& endArgs = loopEnd.getArguments();
        for (size_t i = 0; i < endArgs.size(); i++) {
            assert(endArgs[i].getOperation() != NULL);
            findNonIndexedNodes(*endArgs[i].getOperation(), nonIndexed, loopIndex);
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
                                                            std::set<OperationNode<Base>*>& nonIndexed,
                                                            const Index& loopIndex) {
        if (node.getColor() > 0)
            return node.getColor() == 1;

        if (node.getOperationType() == CGIndexOp) {
            IndexOperationNode<Base>& indexNode = static_cast<IndexOperationNode<Base>&> (node);
            if (&indexNode.getIndex() == &loopIndex) {
                node.setColor(2);
                return false; // depends on the loop index
            }
        }

        const std::vector<Argument<Base> >& args = node.getArguments();
        size_t size = args.size();

        bool indexedPath = false; // whether or not this node depends on indexed independents
        bool nonIndexedArgs = false; // whether or not there are non indexed arguments
        for (size_t a = 0; a < size; a++) {
            OperationNode<Base>* arg = args[a].getOperation();
            if (arg != NULL) {
                bool nonIndexedArg = findNonIndexedNodes(*arg, nonIndexed, loopIndex);
                nonIndexedArgs |= nonIndexedArg;
                indexedPath |= !nonIndexedArg;
            }
        }

        node.setColor(indexedPath ? 2 : 1);

        if (indexedPath && nonIndexedArgs) {
            for (size_t a = 0; a < size; a++) {
                OperationNode<Base>* arg = args[a].getOperation();
                if (arg != NULL && arg->getColor() == 1 && // must be a non indexed expression
                        arg->getOperationType() != CGInvOp) {// no point in moving just one variable outside

                    nonIndexed.insert(arg);
                }
            }
        }

        return !indexedPath;
    }
}

#endif
