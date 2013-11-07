#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_HESS_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_HESS_INCLUDED
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

    namespace loops {

        class HessianElement {
        public:
            size_t location; // location in the compressed hessian vector
            size_t row;
            unsigned short count; // number of times to be added to that location

            inline HessianElement() :
                location(std::numeric_limits<size_t>::max()),
                row(std::numeric_limits<size_t>::max()),
                count(0) {
            }

        };

        template<class Base>
        std::pair<CG<Base>, IndexPattern*> createHessianContribution(CodeHandler<Base>& handler,
                                                                     const std::vector<HessianElement>& positions,
                                                                     const CG<Base>& ddfdxdx,
                                                                     IndexOperationNode<Base>& iterationIndexOp,
                                                                     vector<IfElseInfo<Base> >& ifElses);

    }

    /***************************************************************************
     *  Methods related with loop insertion into the operation graph
     **************************************************************************/

    template<class Base>
    void CLangCompileModelHelper<Base>::analyseSparseHessianWithLoops(const std::vector<size_t>& lowerHessRows,
                                                                      const std::vector<size_t>& lowerHessCols,
                                                                      const std::vector<size_t>& lowerHessOrder,
                                                                      vector<std::set<size_t> >& noLoopEvalJacSparsity,
                                                                      vector<std::set<size_t> >& noLoopEvalHessSparsity,
                                                                      vector<std::map<size_t, std::set<size_t> > >& noLoopEvalHessLocations,
                                                                      std::map<LoopModel<Base>*, loops::HessianWithLoopsInfo<Base> >& loopHessInfo,
                                                                      bool useSymmetry) {
        using namespace std;
        using namespace CppAD::loops;
        using CppAD::vector;

        size_t nonIndexdedEqSize = _funNoLoops != NULL ? _funNoLoops->getOrigDependentIndexes().size() : 0;

        /**
         * determine sparsities
         */
        typename std::set<LoopModel<Base>*>::const_iterator itloop;
        for (itloop = _loopTapes.begin(); itloop != _loopTapes.end(); ++itloop) {
            LoopModel<Base>* l = *itloop;
            l->evalJacobianSparsity();
            l->evalHessianSparsity();
        }

        if (_funNoLoops != NULL) {
            _funNoLoops->evalJacobianSparsity();
            _funNoLoops->evalHessianSparsity();
        }

        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        size_t nnz = lowerHessRows.size();

        noLoopEvalJacSparsity.resize(_funNoLoops != NULL ? m : 0);
        noLoopEvalHessSparsity.resize(_funNoLoops != NULL ? n : 0);
        noLoopEvalHessLocations.resize(noLoopEvalHessSparsity.size());

        loopHessInfo.clear();
        for (itloop = _loopTapes.begin(); itloop != _loopTapes.end(); ++itloop) {
            LoopModel<Base>* loop = *itloop;
            loopHessInfo[loop] = HessianWithLoopsInfo<Base>(*loop);
            loopHessInfo[loop].noLoopEvalHessTempsSparsity.resize(_funNoLoops != NULL ? n : 0);
        }

        /** 
         * Load locations in the compressed hessian
         * d      d y_i
         * d x_j2 d x_j1
         */
        for (size_t eh = 0; eh < nnz; eh++) {
            size_t j1 = lowerHessRows[eh];
            size_t j2 = lowerHessCols[eh];
            size_t e = lowerHessOrder[eh];

            if (_funNoLoops != NULL) {
                // considers only the pattern for the original equations and leaves out the temporaries
                const vector<std::set<size_t> >& dydxx = _funNoLoops->getHessianOrigEqsSparsity();
                if (dydxx.size() > 0) {
                    if (dydxx[j1].find(j2) != dydxx[j1].end()) {
                        /**
                         * Present in the equations outside the loops
                         */
                        noLoopEvalHessSparsity[j1].insert(j2);
                        noLoopEvalHessLocations[j1][j2].insert(e);
                    }
                }
            }

            for (itloop = _loopTapes.begin(); itloop != _loopTapes.end(); ++itloop) {
                LoopModel<Base>* loop = *itloop;
                size_t iterations = loop->getIterationCount();
                const vector<set<size_t> >& loopJac = loop->getJacobianSparsity();
                const vector<set<size_t> >& loopHess = loop->getHessianSparsity();
                HessianWithLoopsInfo<Base>& loopInfo = loopHessInfo.at(loop);

                const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = loop->getIndexedIndepIndexes();
                const std::vector<LoopPosition>& nonIndexedIndepIndexes = loop->getNonIndexedIndepIndexes();
                const std::vector<LoopPosition>& temporaryIndependents = loop->getTemporaryIndependents();

                size_t nIndexed = indexedIndepIndexes.size();
                size_t nNonIndexed = nonIndexedIndepIndexes.size();

                const LoopPosition* posJ1 = loop->getNonIndexedIndepIndexes(j1);
                const LoopPosition* posJ2 = loop->getNonIndexedIndepIndexes(j2);

                /**
                 * indexed - indexed
                 * d      d f_i
                 * d x_l2 d x_l1
                 */
                const std::vector<set<pairss> >& iter2tapeII = loop->getHessianIndexedIndexedTapeIndexes(j1, j2);
                for (size_t iteration = 0; iteration < iter2tapeII.size(); iteration++) {
                    const set<pairss>& tapePairs = iter2tapeII[iteration];

                    set<pairss> ::const_iterator itPairs;
                    for (itPairs = tapePairs.begin(); itPairs != tapePairs.end(); ++itPairs) {
                        size_t tape1 = itPairs->first;
                        size_t tape2 = itPairs->second;
                        pairss tape;
                        if (useSymmetry && tape1 > tape2 && loopHess[tape2].find(tape1) != loopHess[tape2].end()) {
                            tape = pairss(tape2, tape1); // work the symmetry
                        } else {
                            tape = *itPairs;
                        }

                        std::vector<HessianElement>& positions = loopInfo.indexedIndexedPositions[tape];
                        positions.resize(iterations);

                        positions[iteration].location = e;
                        positions[iteration].row = j1;
                        positions[iteration].count++;
                        loopInfo.evalHessSparsity[tape.first].insert(tape.second);
                    }
                }

                /**
                 * indexed - non-indexed 
                 * d      d f_i    ->   d      d f_i
                 * d x_j2 d x_l1        d x_l2 d x_j1
                 */
                if (posJ2 != NULL) {
                    const std::vector<set<size_t> >& iter2tapeJ1OrigJ2 = loop->getHessianIndexedNonIndexedTapeIndexes(j1, j2);
                    for (size_t iteration = 0; iteration < iter2tapeJ1OrigJ2.size(); iteration++) {
                        const set<size_t>& tapeJ1s = iter2tapeJ1OrigJ2[iteration];

                        set<size_t>::const_iterator ittj1;
                        for (ittj1 = tapeJ1s.begin(); ittj1 != tapeJ1s.end(); ++ittj1) {
                            size_t tapeJ1 = *ittj1;

                            std::vector<HessianElement>* positions;
                            if (useSymmetry && loopHess[posJ2->tape].find(tapeJ1) != loopHess[posJ2->tape].end()) {
                                positions = &loopInfo.nonIndexedIndexedPositions[pairss(posJ2->tape, tapeJ1)];
                                loopInfo.evalHessSparsity[posJ2->tape].insert(tapeJ1);
                            } else {
                                positions = &loopInfo.indexedNonIndexedPositions[pairss(tapeJ1, posJ2->tape)];
                                loopInfo.evalHessSparsity[tapeJ1].insert(posJ2->tape);
                            }

                            positions->resize(iterations);
                            (*positions)[iteration].location = e;
                            (*positions)[iteration].row = j1;
                            (*positions)[iteration].count++;
                        }
                    }
                }

                /**
                 * indexed - constant z
                 * d     d f_i    .   d z_k
                 * d z_k d x_l1       d x_j2
                 */
                if (_funNoLoops != NULL) {
                    map<size_t, set<size_t> > iter2tapeJ1 = loop->getIndexedTapeIndexes(j1);
                    map<size_t, set<size_t> >::const_iterator itIter;
                    for (itIter = iter2tapeJ1.begin(); itIter != iter2tapeJ1.end(); ++itIter) {
                        size_t iteration = itIter->first;
                        const set<size_t>& tapeJ1s = itIter->second;

                        set<size_t>::const_iterator itTapeJ1;
                        for (itTapeJ1 = tapeJ1s.begin(); itTapeJ1 != tapeJ1s.end(); ++itTapeJ1) {
                            size_t tapeJ1 = *itTapeJ1;

                            set<size_t>::const_iterator itz = loopHess[tapeJ1].lower_bound(nIndexed + nNonIndexed);

                            pairss pos(tapeJ1, j2);
                            bool used = false;

                            // loop temporary variables
                            for (; itz != loopHess[tapeJ1].end(); ++itz) {
                                size_t tapeJ = *itz;
                                size_t k = temporaryIndependents[tapeJ - nIndexed - nNonIndexed].original;

                                /**
                                 * check if this temporary depends on j2
                                 */
                                const set<size_t>& sparsity = _funNoLoops->getJacobianSparsity()[nonIndexdedEqSize + k];
                                if (sparsity.find(j2) != sparsity.end()) {
                                    noLoopEvalJacSparsity[nonIndexdedEqSize + k].insert(j2); // element required

                                    std::set<size_t>& evals = loopInfo.indexedTempEvals[pos];

                                    used = true;
                                    evals.insert(k);

                                    size_t tapeK = loop->getTempIndepIndexes(k)->tape;
                                    loopInfo.evalHessSparsity[tapeJ1].insert(tapeK);
                                }
                            }

                            if (used) {
                                std::vector<HessianElement>& positions = loopInfo.indexedTempPositions[pos];
                                positions.resize(iterations);

                                positions[iteration].location = e;
                                positions[iteration].row = j1;
                                positions[iteration].count++;
                            }

                        }

                    }

                }

                /**
                 * non-indexed - indexed
                 * d      d f_i
                 * d x_l2 d x_j1
                 */
                if (posJ1 != NULL) {
                    const std::vector<set<size_t> >& iter2TapeJ2 = loop->getHessianNonIndexedIndexedTapeIndexes(j1, j2);
                    for (size_t iteration = 0; iteration < iter2TapeJ2.size(); iteration++) {
                        const set<size_t>& tapeJ2s = iter2TapeJ2[iteration];

                        set<size_t>::const_iterator ittj2;
                        for (ittj2 = tapeJ2s.begin(); ittj2 != tapeJ2s.end(); ++ittj2) {
                            size_t tapeJ2 = *ittj2;

                            std::vector<HessianElement>& positions = loopInfo.nonIndexedIndexedPositions[pairss(posJ1->tape, tapeJ2)];
                            positions.resize(iterations);

                            positions[iteration].location = e;
                            positions[iteration].row = j1;
                            positions[iteration].count++;
                            loopInfo.evalHessSparsity[posJ1->tape].insert(tapeJ2);
                        }
                    }
                }

                /**
                 * non-indexed - non-indexed
                 * d      d f_i
                 * d x_j2 d x_j1
                 */
                bool jInNonIndexed = false;
                pairss orig(j1, j2);

                if (posJ1 != NULL && posJ2 != NULL) {
                    const set<pairss>& orig1orig2 = loop->getHessianNonIndexedNonIndexedIndexes();
                    if (orig1orig2.find(orig) != orig1orig2.end()) {
                        loopInfo.nonIndexedNonIndexedPosition[orig] = e;
                        loopInfo.nonIndexedNonIndexedEvals.insert(orig);

                        loopInfo.evalHessSparsity[posJ1->tape].insert(posJ2->tape);
                        jInNonIndexed = true;
                    }
                }

                /**
                 * non-indexed - temporaries
                 * d     d f_i   .  d z_k
                 * d z_k d x_j1     d x_j2
                 */
                if (_funNoLoops != NULL && posJ1 != NULL) {

                    const set<size_t>& hessRow = loopHess[posJ1->tape];
                    set<size_t>::const_iterator itz = hessRow.lower_bound(nIndexed + nNonIndexed);

                    // loop temporary variables
                    for (; itz != hessRow.end(); ++itz) {
                        size_t tapeJ = *itz;
                        size_t k = temporaryIndependents[tapeJ - nIndexed - nNonIndexed].original;

                        //jacobian of g for k must have j2
                        const set<size_t>& gJacRow = _funNoLoops->getJacobianSparsity()[nonIndexdedEqSize + k];
                        if (gJacRow.find(j2) != gJacRow.end()) {
                            noLoopEvalJacSparsity[nonIndexdedEqSize + k].insert(j2); // element required

                            if (!jInNonIndexed) {
                                CPPADCG_ASSERT_KNOWN(loopInfo.nonIndexedNonIndexedPosition.find(orig) == loopInfo.nonIndexedNonIndexedPosition.end(),
                                                     "Repeated hessian elements requested");

                                loopInfo.nonIndexedNonIndexedPosition[orig] = e;
                                jInNonIndexed = true;
                            }

                            loopInfo.nonIndexedTempEvals[orig].insert(k);

                            size_t tapeK = loop->getTempIndepIndexes(k)->tape;
                            loopInfo.evalHessSparsity[posJ1->tape].insert(tapeK);
                        }

                    }
                }

                /**
                 * temporaries
                 */
                if (_funNoLoops != NULL) {
                    const vector<set<size_t> >& gJac = _funNoLoops->getJacobianSparsity();
                    size_t nk = _funNoLoops->getTemporaryDependentCount();
                    size_t nOrigEq = _funNoLoops->getTapeDependentCount() - nk;

                    const vector<set<size_t> >& dzdxx = _funNoLoops->getHessianTempEqsSparsity();

                    std::set<size_t> usedTapeJ2;

                    for (size_t k1 = 0; k1 < nk; k1++) {
                        if (gJac[nOrigEq + k1].find(j1) == gJac[nOrigEq + k1].end()) {
                            continue;
                        }

                        const LoopPosition* posK1 = loop->getTempIndepIndexes(k1);
                        if (posK1 == NULL) {
                            continue;
                        }

                        /**
                         * temporary - indexed
                         * d     d f_i
                         * d x_l d z_k1
                         */
                        const std::vector<set<size_t> >& iter2TapeJ2 = loop->getHessianTempIndexedTapeIndexes(k1, j2);
                        for (size_t iteration = 0; iteration < iter2TapeJ2.size(); iteration++) {
                            const set<size_t>& tapeJ2s = iter2TapeJ2[iteration];

                            set<size_t>::const_iterator ittj2;
                            for (ittj2 = tapeJ2s.begin(); ittj2 != tapeJ2s.end(); ++ittj2) {
                                size_t tapeJ2 = *ittj2;

                                std::vector<HessianElement>* positions = NULL;

                                if (useSymmetry && loopHess[tapeJ2].find(posK1->tape) != loopHess[tapeJ2].end()) {
                                    if (usedTapeJ2.find(tapeJ2) == usedTapeJ2.end()) {
                                        pairss pos(tapeJ2, j1);
                                        positions = &loopInfo.indexedTempPositions[pos];
                                    }
                                    loopInfo.evalHessSparsity[tapeJ2].insert(posK1->tape);
                                } else {
                                    if (usedTapeJ2.find(tapeJ2) == usedTapeJ2.end()) {
                                        pairss pos(j1, tapeJ2);
                                        positions = &loopInfo.tempIndexedPositions[pos];
                                    }
                                    loopInfo.evalHessSparsity[posK1->tape].insert(tapeJ2);
                                }

                                if (positions != NULL) {
                                    positions->resize(iterations);

                                    (*positions)[iteration].location = e;
                                    (*positions)[iteration].row = j1;
                                    (*positions)[iteration].count++;
                                    usedTapeJ2.insert(tapeJ2);
                                }

                                std::set<size_t>& evals = loopInfo.indexedTempEvals[pairss(tapeJ2, j1)];
                                evals.insert(k1);

                            }
                        }

                        /**
                         * temporary - non-indexed
                         * d      d f_i
                         * d x_j2 d z_k1
                         */
                        if (posJ2 != NULL) {
                            const set<size_t>& hessRow = loop->getHessianSparsity()[posK1->tape];

                            if (hessRow.find(j2) != hessRow.end()) {
                                if (!jInNonIndexed) {
                                    CPPADCG_ASSERT_KNOWN(loopInfo.nonIndexedNonIndexedPosition.find(orig) == loopInfo.nonIndexedNonIndexedPosition.end(),
                                                         "Repeated hessian elements requested");

                                    loopInfo.nonIndexedNonIndexedPosition[orig] = e;
                                    jInNonIndexed = true;
                                }

                                loopInfo.tempNonIndexedEvals[orig].insert(k1);
                                loopInfo.evalHessSparsity[posK1->tape].insert(posJ2->tape);
                            }
                        }
                        /**
                         * temporary - temporary
                         *    d  d f_i    .  d z_k2
                         * d z_k2 d z_k1     d x_j2
                         */
                        // loop hess row
                        const set<size_t>& hessRow = loop->getHessianSparsity()[posK1->tape];
                        set<size_t>::const_iterator itTapeJ2 = hessRow.lower_bound(nIndexed + nNonIndexed);
                        for (; itTapeJ2 != hessRow.end(); ++itTapeJ2) {
                            size_t tapeK2 = *itTapeJ2;
                            size_t k2 = loop->getTemporaryIndependents()[tapeK2 - nIndexed - nNonIndexed].original;

                            const set<size_t>& jacZk2Row = gJac[nOrigEq + k2];
                            if (jacZk2Row.find(j2) != jacZk2Row.end()) { // is this check truly needed?

                                if (!jInNonIndexed) {
                                    CPPADCG_ASSERT_KNOWN(loopInfo.nonIndexedNonIndexedPosition.find(orig) == loopInfo.nonIndexedNonIndexedPosition.end(),
                                                         "Repeated hessian elements requested");

                                    loopInfo.nonIndexedNonIndexedPosition[orig] = e;
                                    jInNonIndexed = true;
                                }

                                loopInfo.tempTempEvals[orig][k1].insert(k2);
                                loopInfo.evalHessSparsity[posK1->tape].insert(tapeK2);
                                noLoopEvalJacSparsity[nOrigEq + k2].insert(j2);
                            }
                        }


                        //
                        noLoopEvalJacSparsity[nOrigEq + k1].insert(j1);

                        /**
                         * temporary - temporary
                         * d f_i   .  d      d z_k1
                         * d z_k1     d x_j2 d x_j1
                         */
                        if (dzdxx[j1].find(j2) != dzdxx[j1].end()) {

                            for (size_t i = 0; i < loopJac.size(); i++) {
                                const set<size_t>& fJacRow = loopJac[i];

                                if (fJacRow.find(posK1->tape) != fJacRow.end()) {

                                    if (!jInNonIndexed) {
                                        CPPADCG_ASSERT_KNOWN(loopInfo.nonIndexedNonIndexedPosition.find(orig) == loopInfo.nonIndexedNonIndexedPosition.end(),
                                                             "Repeated hessian elements requested");

                                        loopInfo.nonIndexedNonIndexedPosition[orig] = e;
                                        jInNonIndexed = true;
                                    }

                                    loopInfo.nonLoopNonIndexedNonIndexed[orig].insert(k1);
                                    loopInfo.evalJacSparsity[i].insert(posK1->tape);
                                    loopInfo.noLoopEvalHessTempsSparsity[j1].insert(j2);
                                }
                            }
                        }
                    }
                }

            }
        }
    }

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::prepareSparseHessianWithLoops(CodeHandler<Base>& handler,
                                                                                   vector<CGBase>& x,
                                                                                   vector<CGBase>& w,
                                                                                   const std::vector<size_t>& lowerHessRows,
                                                                                   const std::vector<size_t>& lowerHessCols,
                                                                                   const std::vector<size_t>& lowerHessOrder,
                                                                                   const std::map<size_t, size_t>& duplicates) {
        using namespace std;
        using namespace CppAD::loops;
        using CppAD::vector;

        handler.setZeroDependents(true);

        size_t nonIndexdedEqSize = _funNoLoops != NULL ? _funNoLoops->getOrigDependentIndexes().size() : 0;

        size_t maxLoc = _hessSparsity.rows.size();
        vector<CGBase> hess(maxLoc);

        vector<set<size_t> > noLoopEvalJacSparsity;
        vector<set<size_t> > noLoopEvalHessSparsity;
        vector<map<size_t, set<size_t> > > noLoopEvalHessLocations;
        map<LoopModel<Base>*, HessianWithLoopsInfo<Base> > loopHessInfo;

        /** 
         * Load locations in the compressed Hessian
         * d      d y_i
         * d x_j2 d x_j1
         */
        analyseSparseHessianWithLoops(lowerHessRows, lowerHessCols, lowerHessOrder,
                                      noLoopEvalJacSparsity, noLoopEvalHessSparsity,
                                      noLoopEvalHessLocations, loopHessInfo, true);

        /***********************************************************************
         *        generate the operation graph
         **********************************************************************/
        /**
         * prepare loop independents
         */
        IndexDclrOperationNode<Base>* iterationIndexDcl = new IndexDclrOperationNode<Base>(LoopModel<Base>::ITERATION_INDEX_NAME);
        handler.manageOperationNodeMemory(iterationIndexDcl);

        // temporaries (zero orders)
        vector<CGBase> tmpsAlias;
        if (_funNoLoops != NULL) {
            ADFun<CGBase>& fun = _funNoLoops->getTape();

            tmpsAlias.resize(fun.Range() - nonIndexdedEqSize);
            for (size_t k = 0; k < tmpsAlias.size(); k++) {
                // to be defined later
                tmpsAlias[k] = handler.createCG(new OperationNode<Base>(CGAliasOp));
            }
        }

        /**
         * prepare loop independents
         */
        typename map<LoopModel<Base>*, HessianWithLoopsInfo<Base> >::iterator itLoop2Info;
        for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info) {
            LoopModel<Base>& lModel = *itLoop2Info->first;
            HessianWithLoopsInfo<Base>& info = itLoop2Info->second;

            /**
             * make the loop start
             */
            info.loopStart = new LoopStartOperationNode<Base>(*iterationIndexDcl, lModel.getIterationCount());
            handler.manageOperationNodeMemory(info.loopStart);

            info.iterationIndexOp = new IndexOperationNode<Base>(*info.loopStart);
            handler.manageOperationNodeMemory(info.iterationIndexOp);
            set<IndexOperationNode<Base>*> indexesOps;
            indexesOps.insert(info.iterationIndexOp);

            /**
             * make the loop's indexed variables
             */
            vector<CGBase> indexedIndeps = createIndexedIndependents(handler, lModel, *info.iterationIndexOp);
            info.x = createLoopIndependentVector(handler, lModel, indexedIndeps, x, tmpsAlias);

            info.w = createLoopDependentVector(handler, lModel, *info.iterationIndexOp);
        }

        /**
         * Calculate Hessians and Jacobians
         */
        /**
         * Loops - evaluate Jacobian and Hessian
         */
        for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info) {
            LoopModel<Base>& lModel = *itLoop2Info->first;
            HessianWithLoopsInfo<Base>& info = itLoop2Info->second;

            _cache.str("");
            _cache << "model (Jacobian + Hessian, loop " << lModel.getLoopId() << ")";
            std::string jobName = _cache.str();
            _cache.str("");
            startingJob("'" + jobName + "'", JobTimer::GRAPH);

            info.evalLoopModelJacobianHessian(false);

            finishedJob();
        }

        /**
         * No loops
         */
        // Jacobian for temporaries
        map<size_t, map<size_t, CGBase> > dzDx;

        if (_funNoLoops != NULL) {
            ADFun<CGBase>& fun = _funNoLoops->getTape();
            vector<CGBase> yNL(fun.Range());

            /**
             * Jacobian and Hessian - temporary variables
             */
            startingJob("'model (Jacobian + Hessian, temporaries)'", JobTimer::GRAPH);

            dzDx = _funNoLoops->calculateJacobianHessianUsedByLoops(loopHessInfo, x, yNL,
                                                                    noLoopEvalJacSparsity,
                                                                    false);

            finishedJob();

            for (size_t i = 0; i < tmpsAlias.size(); i++)
                tmpsAlias[i].getOperationNode()->getArguments().push_back(asArgument(yNL[nonIndexdedEqSize + i]));

            for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info) {
                HessianWithLoopsInfo<Base>& info = itLoop2Info->second;
                // not needed anymore:
                info.dyiDzk.clear();
            }

            /**
             * Hessian - original equations
             */
            _funNoLoops->calculateHessian4OrignalEquations(x, w,
                                                           noLoopEvalHessSparsity, noLoopEvalHessLocations,
                                                           hess);
        }

        /**
         * Loops - Hessian
         */
        for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info) {
            LoopModel<Base>& lModel = *itLoop2Info->first;
            HessianWithLoopsInfo<Base>& info = itLoop2Info->second;

            // store results in indexedLoopResults
            size_t hessElSize = info.indexedIndexedPositions.size() +
                    info.indexedTempPositions.size() +
                    info.nonIndexedIndexedPositions.size() +
                    info.nonIndexedNonIndexedPosition.size();

            if (hessElSize == 0)
                continue; // no second order information

            vector<pair<CGBase, IndexPattern*> > indexedLoopResults(hessElSize);
            size_t hessLE = 0;

            /*******************************************************************
             * indexed - indexed
             */
            map<pairss, std::vector<HessianElement> >::const_iterator it;
            for (it = info.indexedIndexedPositions.begin(); it != info.indexedIndexedPositions.end(); ++it) {
                size_t tapeJ1 = it->first.first;
                size_t tapeJ2 = it->first.second;
                const std::vector<HessianElement>& positions = it->second;

                indexedLoopResults[hessLE++] = createHessianContribution(handler, positions, info.hess[tapeJ1].at(tapeJ2),
                                                                         *info.iterationIndexOp, info.ifElses);
            }

            /**
             * indexed - non-indexed
             * - usually done by  (non-indexed - indexed) by exploiting the symmetry
             */
            for (it = info.indexedNonIndexedPositions.begin(); it != info.indexedNonIndexedPositions.end(); ++it) {
                size_t tapeJ1 = it->first.first;
                size_t tapeJ2 = it->first.second;
                const std::vector<HessianElement>& positions = it->second;

                indexedLoopResults[hessLE++] = createHessianContribution(handler, positions, info.hess[tapeJ1].at(tapeJ2),
                                                                         *info.iterationIndexOp, info.ifElses);
            }

            /**
             * indexed - temporary
             */
            map<pairss, set<size_t> >::const_iterator itEval;
            if (!info.indexedTempPositions.empty()) {
                for (itEval = info.indexedTempEvals.begin(); itEval != info.indexedTempEvals.end(); ++itEval) {
                    size_t tapeJ1 = itEval->first.first;
                    size_t j2 = itEval->first.second;
                    const set<size_t>& ks = itEval->second;

                    map<pairss, std::vector<HessianElement> >::const_iterator itPos = info.indexedTempPositions.find(itEval->first);
                    if (itPos != info.indexedTempPositions.end()) {
                        const std::vector<HessianElement>& positions = itPos->second;

                        CGBase hessVal = Base(0);
                        set<size_t>::const_iterator itz;
                        for (itz = ks.begin(); itz != ks.end(); ++itz) {
                            size_t k = *itz;
                            size_t tapeK = lModel.getTempIndepIndexes(k)->tape;
                            hessVal += info.hess[tapeJ1].at(tapeK) * dzDx[k][j2];
                        }

                        indexedLoopResults[hessLE++] = createHessianContribution(handler, positions, hessVal,
                                                                                 *info.iterationIndexOp, info.ifElses);
                    }
                }
            }

            /*******************************************************************
             * non-indexed - indexed
             */
            for (it = info.nonIndexedIndexedPositions.begin(); it != info.nonIndexedIndexedPositions.end(); ++it) {
                size_t tapeJ1 = it->first.first;
                size_t tapeJ2 = it->first.second;
                const std::vector<HessianElement>& positions = it->second;

                indexedLoopResults[hessLE++] = createHessianContribution(handler, positions, info.hess[tapeJ1].at(tapeJ2),
                                                                         *info.iterationIndexOp, info.ifElses);
            }

            /*******************************************************************
             * temporary - indexed
             * 
             *      d f_i       .    d x_k1
             * d x_l2  d z_k1        d x_j1
             * 
             * -> usually done by  (indexed - temporary) by exploiting the symmetry
             */
            if (!info.tempIndexedPositions.empty()) {
                for (itEval = info.indexedTempEvals.begin(); itEval != info.indexedTempEvals.end(); ++itEval) {
                    size_t tapeJ2 = itEval->first.first;
                    size_t j1 = itEval->first.second;
                    const set<size_t>& ks = itEval->second;

                    map<pairss, std::vector<HessianElement> >::const_iterator itPos = info.tempIndexedPositions.find(pairss(j1, tapeJ2));
                    if (itPos != info.tempIndexedPositions.end()) {
                        const std::vector<HessianElement>& positions = itPos->second;
                        CGBase hessVal = Base(0);
                        set<size_t>::const_iterator itz;
                        for (itz = ks.begin(); itz != ks.end(); ++itz) {
                            size_t k = *itz;
                            size_t tapeK = lModel.getTempIndepIndexes(k)->tape;
                            hessVal += info.hess[tapeK].at(tapeJ2) * dzDx[k][j1];
                        }

                        indexedLoopResults[hessLE++] = createHessianContribution(handler, positions, hessVal,
                                                                                 *info.iterationIndexOp, info.ifElses);
                    }
                }
            }

            /*******************************************************************
             * contributions to a constant location
             */
            map<pairss, size_t>::const_iterator orig2PosIt;
            for (orig2PosIt = info.nonIndexedNonIndexedPosition.begin(); orig2PosIt != info.nonIndexedNonIndexedPosition.end(); ++orig2PosIt) {
                const pairss& orig = orig2PosIt->first;
                size_t e = orig2PosIt->second;

                size_t j1 = orig.first;
                size_t j2 = orig.second;
                const LoopPosition* posJ1 = lModel.getNonIndexedIndepIndexes(j1);
                const LoopPosition* posJ2 = lModel.getNonIndexedIndepIndexes(j2);

                // location
                IndexPattern* pattern = new LinearIndexPattern(0, 0, 0, e);
                handler.manageLoopDependentIndexPattern(pattern);

                /**
                 * non-indexed - non-indexed
                 */
                CGBase hessVal = Base(0);

                if (info.nonIndexedNonIndexedEvals.find(orig) != info.nonIndexedNonIndexedEvals.end()) {
                    hessVal = info.hess[posJ1->tape].at(posJ2->tape);
                }

                /**
                 * non-indexed - temporary
                 */
                map<pairss, set<size_t> >::const_iterator itNT = info.nonIndexedTempEvals.find(orig);
                if (itNT != info.nonIndexedTempEvals.end()) {
                    const set<size_t>& ks = itNT->second;

                    set<size_t>::const_iterator itz;
                    for (itz = ks.begin(); itz != ks.end(); ++itz) {
                        size_t k = *itz;
                        size_t tapeK = lModel.getTempIndepIndexes(k)->tape;
                        hessVal += info.hess[posJ1->tape].at(tapeK) * dzDx[k][j2];
                    }
                }

                /**
                 * temporary - non-indexed 
                 * 
                 *      d f_i       .    d x_k1
                 * d x_j2  d z_k1        d x_j1
                 */
                map<pairss, set<size_t> >::const_iterator itTN = info.tempNonIndexedEvals.find(orig);
                if (itTN != info.tempNonIndexedEvals.end()) {
                    const set<size_t>& ks = itTN->second;

                    set<size_t>::const_iterator itz;
                    for (itz = ks.begin(); itz != ks.end(); ++itz) {
                        size_t k1 = *itz;
                        size_t tapeK = lModel.getTempIndepIndexes(k1)->tape;
                        hessVal += info.hess[tapeK].at(posJ2->tape) * dzDx[k1][j1];
                    }
                }

                /**
                 * temporary - temporary
                 */
                map<pairss, map<size_t, set<size_t> > >::const_iterator itTT = info.tempTempEvals.find(orig);
                if (itTT != info.tempTempEvals.end()) {
                    const map<size_t, set<size_t> >& k1k2 = itTT->second;

                    CGBase sum = Base(0);

                    map<size_t, set<size_t> >::const_iterator itzz;
                    for (itzz = k1k2.begin(); itzz != k1k2.end(); ++itzz) {
                        size_t k1 = itzz->first;
                        const set<size_t>& k2s = itzz->second;
                        size_t tapeK1 = lModel.getTempIndepIndexes(k1)->tape;

                        CGBase tmp = Base(0);
                        for (set<size_t>::const_iterator itk2 = k2s.begin(); itk2 != k2s.end(); ++itk2) {
                            size_t k2 = *itk2;
                            size_t tapeK2 = lModel.getTempIndepIndexes(k2)->tape;

                            tmp += info.hess[tapeK1].at(tapeK2) * dzDx[k2][j2];
                        }

                        sum += tmp * dzDx[k1][j1];
                    }

                    hessVal += sum;
                }

                /**
                 * temporary - temporary
                 */
                map<pairss, set<size_t> >::const_iterator itTT2 = info.nonLoopNonIndexedNonIndexed.find(orig);
                if (itTT2 != info.nonLoopNonIndexedNonIndexed.end()) {
                    hessVal += info.dzDxx.at(j1).at(j2); // it is already the sum of ddz / dx_j1 dx_j2
                }

                // place the result
                indexedLoopResults[hessLE++] = make_pair(hessVal, pattern);
            }


            assert(hessLE == indexedLoopResults.size());

            /**
             * make the loop end
             */
            size_t assignOrAdd = 1;
            set<IndexOperationNode<Base>*> indexesOps;
            indexesOps.insert(info.iterationIndexOp);
            info.loopEnd = createLoopEnd(handler, *info.loopStart, indexedLoopResults, indexesOps, assignOrAdd);

            std::vector<size_t>::const_iterator itE;
            for (itE = lowerHessOrder.begin(); itE != lowerHessOrder.end(); ++itE) {
                // an additional alias variable is required so that each dependent variable can have its own ID
                size_t e = *itE;
                if (hess[e].isParameter() && hess[e].IdenticalZero()) {
                    hess[e] = handler.createCG(new OperationNode<Base> (CGDependentMultiAssignOp, Argument<Base>(*info.loopEnd)));

                } else if (hess[e].getOperationNode() != NULL && hess[e].getOperationNode()->getOperationType() == CGDependentMultiAssignOp) {
                    hess[e].getOperationNode()->getArguments().push_back(Argument<Base>(*info.loopEnd));

                } else {
                    hess[e] = handler.createCG(new OperationNode<Base> (CGDependentMultiAssignOp, asArgument(hess[e]), Argument<Base>(*info.loopEnd)));
                }
            }

            // not needed anymore:
            info.hess.clear();
            info.dzDxx.clear();

            /**
             * move non-indexed expressions outside loop
             */
            moveNonIndexedOutsideLoop(*info.loopStart, *info.loopEnd);
        }

        /**
         * duplicates (TODO: use loops)
         */
        // make use of the symmetry of the Hessian in order to reduce operations
        std::map<size_t, size_t>::const_iterator it2;
        for (it2 = duplicates.begin(); it2 != duplicates.end(); ++it2) {
            if (hess[it2->second].isVariable())
                hess[it2->first] = handler.createCG(new OperationNode<Base> (CGAliasOp, asArgument(hess[it2->second])));
            else
                hess[it2->first] = hess[it2->second].getValue();
        }

        return hess;
    }

    namespace loops {

        template<class Base>
        std::pair<CG<Base>, IndexPattern*> createHessianContribution(CodeHandler<Base>& handler,
                                                                     const std::vector<HessianElement>& positions,
                                                                     const CG<Base>& ddfdxdx,
                                                                     IndexOperationNode<Base>& iterationIndexOp,
                                                                     vector<IfElseInfo<Base> >& ifElses) {
            using namespace std;

            // combine iterations with the same number of additions
            map<size_t, map<size_t, size_t> > locations;
            for (size_t iter = 0; iter < positions.size(); iter++) {
                size_t c = positions[iter].count;
                if (c > 0) {
                    locations[c][iter] = positions[iter].location;
                }
            }

            map<size_t, CG<Base> > results;

            // generate the index pattern for the hessian compressed element
            map<size_t, map<size_t, size_t> >::const_iterator countIt;
            for (countIt = locations.begin(); countIt != locations.end(); ++countIt) {
                size_t count = countIt->first;

                CG<Base> val = ddfdxdx;
                for (size_t c = 1; c < count; c++)
                    val += ddfdxdx;

                results[count] = val;
            }

            if (results.size() == 1 && locations.begin()->second.size() == positions.size()) {
                // same expression present in all iterations

                // generate the index pattern for the Hessian compressed element
                IndexPattern* pattern = IndexPattern::detect(locations.begin()->second);
                handler.manageLoopDependentIndexPattern(pattern);

                return make_pair(results.begin()->second, pattern);

            } else {
                /**
                 * must create a conditional element so that this 
                 * contribution to the hessian is only evaluated at the
                 * relevant iterations
                 */
                map<size_t, map<size_t, size_t> >::const_iterator countIt;

                // try to find an existing if-else where these operations can be added
                map<SizeN1stIt, pair<size_t, set<size_t> > > firstIt2Count2Iterations;
                for (countIt = locations.begin(); countIt != locations.end(); ++countIt) {
                    set<size_t> iterations;
                    mapKeys(countIt->second, iterations);
                    SizeN1stIt pos(iterations.size(), *iterations.begin());
                    firstIt2Count2Iterations[pos] = make_pair(countIt->first, iterations);
                }

                IfElseInfo<Base>* ifElseBranches = findExistingIfElse(ifElses, firstIt2Count2Iterations);
                bool reusingIfElse = ifElseBranches != NULL;
                if (!reusingIfElse) {
                    size_t s = ifElses.size();
                    ifElses.resize(s + 1);
                    ifElseBranches = &ifElses[s];
                }

                /**
                 * create/change each if/else branch
                 */
                OperationNode<Base>* ifStart = NULL;
                OperationNode<Base>* ifBranch = NULL;
                Argument<Base> nextBranchArg;
                set<size_t> usedIter;

                map<SizeN1stIt, pair<size_t, set<size_t> > >::const_iterator it1st2Count2Iters;
                for (it1st2Count2Iters = firstIt2Count2Iterations.begin(); it1st2Count2Iters != firstIt2Count2Iterations.end(); ++it1st2Count2Iters) {
                    size_t firstIt = it1st2Count2Iters->first.second;
                    size_t count = it1st2Count2Iters->second.first;
                    const set<size_t>& iterations = it1st2Count2Iters->second.second;

                    size_t iterCount = iterations.size();

                    SizeN1stIt pos(iterCount, firstIt);

                    if (reusingIfElse) {
                        //reuse existing node
                        ifBranch = ifElseBranches->firstIt2Branch.at(pos).node;
                        if (nextBranchArg.getOperation() != NULL)
                            ifBranch->getArguments().push_back(nextBranchArg);

                    } else if (usedIter.size() + iterCount == positions.size()) {
                        // all other iterations: ELSE
                        ifBranch = new OperationNode<Base>(CGElseOp, Argument<Base>(*ifBranch), nextBranchArg);
                        handler.manageOperationNodeMemory(ifBranch);
                    } else {
                        // depends on the iterations indexes
                        OperationNode<Base>* cond = createIndexConditionExpressionOp<Base>(iterations, usedIter, positions.size() - 1, iterationIndexOp);
                        handler.manageOperationNodeMemory(cond);

                        if (ifStart == NULL) {
                            // IF
                            ifStart = new OperationNode<Base>(CGStartIfOp, Argument<Base>(*cond));
                            ifBranch = ifStart;
                        } else {
                            // ELSE IF
                            ifBranch = new OperationNode<Base>(CGElseIfOp, Argument<Base>(*ifBranch), Argument<Base>(*cond), nextBranchArg);
                        }

                        handler.manageOperationNodeMemory(ifBranch);

                        usedIter.insert(iterations.begin(), iterations.end());
                    }

                    const map<size_t, size_t>& locationsC = locations[count];
                    IndexPattern* pattern = IndexPattern::detect(locationsC);
                    handler.manageLoopDependentIndexPattern(pattern);

                    std::vector<size_t> ainfo(2);
                    ainfo[0] = handler.addLoopDependentIndexPattern(*pattern); // dependent index pattern location
                    ainfo[1] = 1; // assignOrAdd
                    std::vector<Argument<Base> > indexedArgs(2);
                    indexedArgs[0] = asArgument(results[count]); // indexed expression
                    indexedArgs[1] = Argument<Base>(iterationIndexOp); // dependency on the index
                    OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, ainfo, indexedArgs);
                    handler.manageOperationNodeMemory(yIndexed);

                    OperationNode<Base>* ifAssign = new OperationNode<Base>(CGCondResultOp, Argument<Base>(*ifBranch), Argument<Base>(*yIndexed));
                    handler.manageOperationNodeMemory(ifAssign);
                    nextBranchArg = Argument<Base>(*ifAssign);

                    if (!reusingIfElse) {
                        IfBranchInfo<Base>& branch = ifElseBranches->firstIt2Branch[pos]; // creates a new if branch
                        branch.iterations = iterations;
                        branch.node = ifBranch;
                    }
                }

                /**
                 * end if
                 */
                if (reusingIfElse) {
                    ifElseBranches->endIf->getArguments().push_back(nextBranchArg);
                } else {
                    ifElseBranches->endIf = new OperationNode<Base>(CGEndIfOp, Argument<Base>(*ifBranch), nextBranchArg);
                    handler.manageOperationNodeMemory(ifElseBranches->endIf);
                }

                IndexPattern* p = NULL;
                return make_pair(handler.createCG(Argument<Base>(*ifElseBranches->endIf)), p);
            }

        }

    }
}

#endif
