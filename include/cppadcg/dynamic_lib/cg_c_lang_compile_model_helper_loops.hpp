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

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareSparseJacobianWithLoops(CodeHandler<Base>& handler,
                                                                       std::vector<CGBase>& jac) {
        size_t nnz = _jacSparsity.rows.size();
        // loop -> equation pattern -> tape independent -> orig independent(temporaries only) -> iteration = position
        std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > > jacIndexPatterns;

        std::vector<JacOrigElementLoopInfo<Base> > garbageCollection;
        garbageCollection.reserve(nnz);

        std::map<OperationNode<Base>*, OperationNode<Base>*> loopDependents;
        
        for (size_t e = 0; e < nnz; e++) {
            OperationNode<Base>* jacNode = jac[e].getOperationNode();
            if (jacNode != NULL && jacNode->getOperationType() == CGLoopResultOp) {
                OperationNode<Base>* loop = handleSparseJacLoopResult(handler, e, nnz, jacNode, jacIndexPatterns, Base(1.0), garbageCollection);
                loopDependents[jacNode] = loop;
            } else {
                OperationNode<Base>* loop = findSparseJacLoopResult(handler, e, nnz, jacNode, jacIndexPatterns, garbageCollection);
                if (loop != NULL) {
                    loopDependents[jacNode] = loop;
                }
            }
        }

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
                        JacOrigElementLoopInfo<Base>* orig = tapeEleInfo.origJacElement;

                        // make sure all element are requested
                        std::vector<size_t>::const_iterator ite;
                        for (ite = orig->jacIndexes.begin(); ite != orig->jacIndexes.end(); ++ite) {
                            if (*ite == nnz) {
                                throw CGException("All jacobian elements of an equation pattern (equation in a loop) must be requested for all iterations");
                            }
                        }

                        orig->jacPattern = IndexPattern::detect(orig->jacIndexes);
                        // TODO: consider clearing origEleInfo.jacIndexes
                    }
                }
            }
        }

        handler.prepareLoops(loopDependents, jacIndexPatterns);
    }

    template<class Base>
    OperationNode<Base>* CLangCompileModelHelper<Base>::findSparseJacLoopResult(CodeHandler<Base>& handler,
                                                                                size_t e,
                                                                                size_t nnz,
                                                                                OperationNode<Base>* jacNode,
                                                                                std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > >& jacIndexPatterns,
                                                                                std::vector<JacOrigElementLoopInfo<Base> >& garbageCollection) {
        // this is only used for the reverse mode
        if (jacNode == NULL)
            return NULL;

        const std::vector<Argument<Base> >& args = jacNode->getArguments();
        if (jacNode->getOperationType() == CGAddOp) {
            OperationNode<Base>* loop1 = findSparseJacLoopResult(handler, e, nnz, args[0].getOperation(), jacIndexPatterns, garbageCollection);
            OperationNode<Base>* loop2 = findSparseJacLoopResult(handler, e, nnz, args[1].getOperation(), jacIndexPatterns, garbageCollection);
            assert(loop1 == NULL || loop2 == NULL || loop1 == loop2);
            if (loop1 != NULL)
                return loop1;
            else
                return loop2;
        } else if (jacNode->getOperationType() == CGMulOp) {
            OperationNode<Base>* resultNode = args[0].getOperation(); // TODO: maybe also consider args[1]
            if (resultNode != NULL && resultNode->getOperationType() == CGLoopResultOp) {
                return handleSparseJacLoopResult(handler, e, nnz, resultNode, jacIndexPatterns, args[1], garbageCollection);
            }
        } else if (jacNode->getOperationType() == CGLoopResultOp) {
            return handleSparseJacLoopResult(handler, e, nnz, jacNode, jacIndexPatterns, Base(1.0), garbageCollection);
        }

        return NULL;
    }

    template<class Base>
    OperationNode<Base>* CLangCompileModelHelper<Base>::handleSparseJacLoopResult(CodeHandler<Base>& handler,
                                                                                  size_t e,
                                                                                  size_t nnz,
                                                                                  OperationNode<Base>* node,
                                                                                  std::map<LoopAtomicFun<Base>*, std::map<size_t, std::map<size_t, JacTapeElementLoopInfo<Base> > > >& jacIndexPatterns,
                                                                                  Argument<Base> arg,
                                                                                  std::vector<JacOrigElementLoopInfo<Base> >& garbageCollection) {

        OperationNode<Base>* loopNode = node->getArguments()[0].getOperation();
        CGOpCode loopOpType = loopNode->getOperationType();
        size_t loopId = loopNode->getInfo()[0];
        //size_t p = loopNode->getInfo()[3];
        assert((loopOpType == CGLoopForwardOp && loopNode->getInfo()[3] == 1) ||
               (loopOpType == CGLoopReverseOp && loopNode->getInfo()[3] == 0));

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

        JacOrigElementLoopInfo<Base>* origJacEl = NULL;
        JacTapeElementLoopInfo<Base>& ref = lRawIndexes[tapeI][*tapeJs.begin()];
        if (ref.origIndep2Info.size() != 0) {
            bool isTemporary = loop->isTemporary(*tapeJs.begin());
            size_t jRef = 0;
            if (isTemporary) {
                jRef = j;
            }
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
        }
        origJacEl->jacIndexes[iteration] = e;


        std::set<size_t>::const_iterator itTapeJ;
        for (itTapeJ = tapeJs.begin(); itTapeJ != tapeJs.end(); ++itTapeJ) {
            size_t tapeJ = *itTapeJ;

            JacTapeElementLoopInfo<Base>& rawLoopIJ = lRawIndexes[tapeI][tapeJ];
            JacOrigElementIndepLoopInfo<Base>* tapeJacElement;
            // check if it is a temporary
            if (loop->isTemporary(tapeJ)) {
                size_t nIndexed = loop->getIndexedIndepIndexes().size();
                size_t nNonIndexed = loop->getNonIndexedIndepIndexes().size();
                const std::vector<LoopPositionTmp>& tmps = loop->getTemporaryIndependents();
                const LoopPositionTmp& pos = tmps[tapeJ - nIndexed - nNonIndexed];

                tapeJacElement = &rawLoopIJ.origIndep2Info[j];
                if (loopOpType == CGLoopForwardOp)
                    tapeJacElement->arg = loopNode->getArguments()[pos.atomic * 2 + 1]; // args = [tx]
                else //CGLoopReverseOp
                    tapeJacElement->arg = arg; // args = [ x , py ] ////////////////////// TODO!!!

            } else {
                tapeJacElement = &rawLoopIJ.origIndep2Info[0];
                tapeJacElement->arg = Argument<Base>(Base(1));
            }

            tapeJacElement->origJacElement = origJacEl;
        }

        return loopNode;
    }

}

#endif
