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

    namespace loops {

        class HessianElement {
        public:
            size_t location; // location in the compressed hessian vector
            unsigned short count; // number of times to be added to that location

            inline HessianElement() :
                location(std::numeric_limits<size_t>::max()),
                count(0) {
            }

        };
        
        template<class Base>
        class HessianWithLoopsInfo;
        
        template<class Base>
        std::pair<CG<Base>, IndexPattern*> createHessianContribution(CodeHandler<Base>& handler,
                                                                     const std::vector<HessianElement>& positions,
                                                                     const CG<Base>& ddfdxdx,
                                                                     IndexOperationNode<Base>& iterationIndexOp);

        template<class Base>
        OperationNode<Base>* createIndexConditionExpression(const std::set<size_t>& iterations,
                                                            const std::set<size_t>& usedIter,
                                                            size_t maxIter,
                                                            IndexOperationNode<Base>& iterationIndexOp);

    }

    /***************************************************************************
     *  Methods related with loop insertion into the operation graph
     **************************************************************************/

    template<class Base>
    vector<CG<Base> > CLangCompileModelHelper<Base>::prepareSparseHessianWithLoops(CodeHandler<Base>& handler,
                                                                                   vector<CGBase>& x,
                                                                                   vector<CGBase>& w,
                                                                                   const std::vector<size_t>& lowerHessRows,
                                                                                   const std::vector<size_t>& lowerHessCols,
                                                                                   const std::vector<size_t>& lowerHessOrder,
                                                                                   const std::map<size_t, size_t>& duplicates) {
        typedef std::pair<size_t, size_t> pairss;

        using namespace std;
        using namespace CppAD::loops;
        using CppAD::vector;

        handler.setZeroDependents(true);

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

        size_t maxLoc = _hessSparsity.rows.size();
        size_t nnz = lowerHessRows.size();
        vector<CGBase> hess(maxLoc);

        vector<set<size_t> > noLoopEvalJacSparsity(_funNoLoops != NULL ? m : 0);

        vector<set<size_t> > noLoopEvalHessSparsity(_funNoLoops != NULL ? n : 0);
        vector<map<size_t, set<size_t> > > noLoopEvalHessLocations(noLoopEvalHessSparsity.size());

        vector<set<size_t> > noLoopEvalHessTempsSparsity(_funNoLoops != NULL ? n : 0);

        map<LoopModel<Base>*, HessianWithLoopsInfo<Base> > loopHessInfo;
        for (itloop = _loopTapes.begin(); itloop != _loopTapes.end(); ++itloop) {
            LoopModel<Base>* loop = *itloop;
            loopHessInfo[loop] = HessianWithLoopsInfo<Base>(*loop);
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
                const set<size_t>& row = _funNoLoops->getHessianSparsity()[j1];
                if (row.find(j2) != row.end()) {
                    /**
                     * Present in the equations outside the loops
                     */
                    noLoopEvalHessSparsity[j1].insert(j2);
                    noLoopEvalHessLocations[j1][j2].insert(e);
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
                        pairss tape; // work the symmetry
                        if (itPairs->first < itPairs->second) {
                            tape = *itPairs;
                        } else {
                            tape = pairss(itPairs->second, itPairs->first);
                        }

                        std::vector<HessianElement>& positions = loopInfo.indexedIndexedPositions[tape];
                        positions.resize(iterations);
                        //CPPADCG_ASSERT_KNOWN(positions[iteration] == maxLoc, "Repeated hessian elements requested");

                        positions[iteration].location = e;
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

                            std::vector<HessianElement>& positions = loopInfo.nonIndexedIndexedPositions[pairss(posJ2->tape, tapeJ1)];
                            positions.resize(iterations);

                            positions[iteration].location = e;
                            positions[iteration].count++;
                            loopInfo.evalHessSparsity[posJ2->tape].insert(tapeJ1);
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

                            // loop temporary variables
                            for (; itz != loopHess[tapeJ1].end(); ++itz) {
                                size_t k = temporaryIndependents[*itz - nIndexed - nNonIndexed].original;

                                /**
                                 * check if this temporary depends on j2
                                 */
                                const set<size_t>& sparsity = _funNoLoops->getJacobianSparsity()[nonIndexdedEqSize + k];
                                if (sparsity.find(j2) != sparsity.end()) {
                                    noLoopEvalJacSparsity[nonIndexdedEqSize + k].insert(j2); // element required

                                    pairss pos(tapeJ1, j2);

                                    std::set<size_t>& evals = loopInfo.indexedTempEvals[pos];
                                    std::vector<HessianElement>& positions = loopInfo.indexedTempPositions[pos];
                                    positions.resize(iterations);

                                    positions[iteration].location = e;
                                    positions[iteration].count++;
                                    evals.insert(k);

                                    size_t tapeK = loop->getTempIndepIndexes(k)->tape;
                                    loopInfo.evalHessSparsity[tapeJ1].insert(tapeK);
                                }
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
                        size_t k = temporaryIndependents[*itz - nIndexed - nNonIndexed].original;

                        //jacobian of g for k must have j2
                        const set<size_t>& gJacRow = loopJac[k];
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

                                pairss pos(tapeJ2, j1);

                                std::set<size_t>& evals = loopInfo.indexedTempEvals[pos];
                                std::vector<HessianElement>& positions = loopInfo.indexedTempPositions[pos];
                                positions.resize(iterations);

                                positions[iteration].location = e;
                                positions[iteration].count++;
                                evals.insert(k1);

                                loopInfo.evalHessSparsity[tapeJ2].insert(posK1->tape);
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
                        const set<size_t>& gHessRow = _funNoLoops->getHessianSparsity()[j1];
                        if (gHessRow.find(j2) != gHessRow.end()) {

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

                                    noLoopEvalHessTempsSparsity[j1].insert(j2);
                                }
                            }
                        }
                    }
                }

            }
        }

        /***********************************************************************
         *        generate the operation graph
         **********************************************************************/
        /**
         * Calculate hessians and jacobians
         */
        // temporaries (zero orders)
        vector<CGBase> tmps;

        /**
         * No loops - zero order
         */
        if (_funNoLoops != NULL) {
            ADFun<CGBase>& fun = _funNoLoops->getTape();
            vector<CGBase> depNL = fun.Forward(0, x);

            tmps.resize(depNL.size() - nonIndexdedEqSize);
            for (size_t i = 0; i < tmps.size(); i++)
                tmps[i] = depNL[nonIndexdedEqSize + i];
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
            info.loopStart = new LoopStartOperationNode<Base>(lModel);
            handler.manageOperationNodeMemory(info.loopStart);

            info.iterationIndexOp = new IndexOperationNode<Base>(LoopModel<Base>::ITERATION_INDEX, *info.loopStart);
            handler.manageOperationNodeMemory(info.iterationIndexOp);
            set<IndexOperationNode<Base>*> indexesOps;
            indexesOps.insert(info.iterationIndexOp);

            vector<CGBase> indexedIndeps = createIndexedIndependents(handler, lModel, *info.iterationIndexOp);
            info.x = createLoopIndependentVector(handler, lModel, indexedIndeps, x, tmps);

            info.w = createLoopDependentVector(handler, lModel, *info.iterationIndexOp);
        }

        /**
         * No loops - Jacobian
         */
        // jacobian for equations outside loops
        vector<CGBase> jacNoLoop;
        // jacobian for temporaries
        map<size_t, map<size_t, CGBase> > dzDx;

        if (_funNoLoops != NULL) {
            ADFun<CGBase>& fun = _funNoLoops->getTape();

            std::vector<size_t> row, col;
            generateSparsityIndexes(noLoopEvalJacSparsity, row, col);

            if (row.size() > 0) {
                jacNoLoop.resize(row.size());

                CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
                if (estimateBestJacobianADMode(row, col)) {
                    fun.SparseJacobianForward(x, _funNoLoops->getJacobianSparsity(), row, col, jacNoLoop, work);
                } else {
                    fun.SparseJacobianReverse(x, _funNoLoops->getJacobianSparsity(), row, col, jacNoLoop, work);
                }

                for (size_t el = 0; el < row.size(); el++) {
                    size_t inl = row[el];
                    size_t j = col[el];
                    assert(inl >= nonIndexdedEqSize);

                    // dz_k/dx_v (for temporary variable)
                    size_t k = inl - nonIndexdedEqSize;
                    dzDx[k][j] = jacNoLoop[el];
                }
            }
        }

        /**
         * Loops - Jacobian
         */
        for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info) {
            LoopModel<Base>& lModel = *itLoop2Info->first;
            HessianWithLoopsInfo<Base>& info = itLoop2Info->second;
            ADFun<CGBase>& fun = lModel.getTape();
            const vector<std::set<size_t> >& jacTapeSparsity = lModel.getJacobianSparsity();

            printSparsityPattern(jacTapeSparsity, "jac - loop");
            printSparsityPattern(info.evalJacSparsity, "jac - loop -eval");

            /**
             * evaluate loop model jacobian
             */
            std::vector<size_t> row, col;
            generateSparsityIndexes(info.evalJacSparsity, row, col);
            if (row.size() > 0) {
                vector<CGBase> jacLoop(row.size());

                CppAD::sparse_jacobian_work work; // temporary structure for CppAD
                if (estimateBestJacobianADMode(row, col)) {
                    fun.SparseJacobianForward(info.x, jacTapeSparsity, row, col, jacLoop, work);
                } else {
                    fun.SparseJacobianReverse(info.x, jacTapeSparsity, row, col, jacLoop, work);
                }

                // save/organize results
                for (size_t el = 0; el < jacLoop.size(); el++) {
                    size_t tapeI = row[el];
                    size_t tapeJ = col[el];
                    info.dyiDzk[tapeI][tapeJ] = jacLoop[el];
                }
            }
        }

        /**
         * No Loops - Hessian
         */
        map<size_t, map<size_t, CGBase> > dzDxx;
        if (_funNoLoops != NULL) {
            ADFun<CGBase>& fun = _funNoLoops->getTape();

            vector<CGBase> wNoLoop(_funNoLoops->getTapeDependentCount());

            vector<CGBase> hessNoLoop;

            const std::vector<size_t>& origIndexes = _funNoLoops->getOrigDependentIndexes();

            /**
             * hessian - original equations
             */
            std::vector<size_t> row, col;
            generateSparsityIndexes(noLoopEvalHessSparsity, row, col);

            if (row.size() > 0) {
                hessNoLoop.resize(row.size());

                for (size_t inl = 0; inl < origIndexes.size(); inl++) {
                    wNoLoop[inl] = w[origIndexes[inl]];
                }

                CppAD::sparse_hessian_work work; // temporary structure for CPPAD
                fun.SparseHessian(x, wNoLoop, _funNoLoops->getHessianSparsity(), row, col, hessNoLoop, work);

                // save non-indexed hessian elements
                for (size_t el = 0; el < row.size(); el++) {
                    size_t j1 = row[el];
                    size_t j2 = col[el];
                    const set<size_t>& locations = noLoopEvalHessLocations[j1][j2];
                    for (set<size_t>::const_iterator itE = locations.begin(); itE != locations.end(); ++itE)
                        hess[*itE] = hessNoLoop[el];
                }
            }

            /**
             * hessian - temporary variables
             */
            generateSparsityIndexes(noLoopEvalHessTempsSparsity, row, col);

            if (row.size() > 0) {
                hessNoLoop.resize(row.size());

                for (size_t inl = 0; inl < origIndexes.size(); inl++) {
                    wNoLoop[inl] = Base(0);
                }

                for (size_t inl = origIndexes.size(); inl < wNoLoop.size(); inl++) {
                    CGBase sum = Base(0);
                    for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info) {
                        HessianWithLoopsInfo<Base>& info = itLoop2Info->second;

                        for (size_t i = 0; i < info.dyiDzk.size(); i++) {
                            const map<size_t, CGBase>& row = info.dyiDzk[i];
                            typename map<size_t, CGBase>::const_iterator itCol;
                            for (itCol = row.begin(); itCol != row.end(); ++itCol) {
                                const CGBase& val = itCol->second;
                                sum += val * info.w[i];
                            }
                        }
                    }

                    wNoLoop[inl] = sum;
                }

                CppAD::sparse_hessian_work workTemps;
                fun.SparseHessian(x, wNoLoop, _funNoLoops->getHessianSparsity(), row, col, hessNoLoop, workTemps);

                // save hessian
                for (size_t el = 0; el < row.size(); el++) {
                    size_t j1 = row[el];
                    size_t j2 = col[el];
                    dzDxx[j1][j2] = hessNoLoop[el];
                }
            }
        }

        /**
         * Loops - Hessian
         */
        for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info) {
            LoopModel<Base>& lModel = *itLoop2Info->first;
            HessianWithLoopsInfo<Base>& info = itLoop2Info->second;
            ADFun<CGBase>& fun = lModel.getTape();

            std::vector<size_t> row, col;
            generateSparsityIndexes(info.evalHessSparsity, row, col);

            if (row.empty())
                continue;

            vector<CGBase> hessLoopFlat(row.size());

            CppAD::sparse_hessian_work work; // temporary structure for CPPAD
            fun.SparseHessian(info.x, info.w, lModel.getHessianSparsity(), row, col, hessLoopFlat, work);

            // save non-indexed hessian elements
            vector<map<size_t, CGBase> > hessLoop(fun.Domain());
            for (size_t el = 0; el < row.size(); el++) {
                size_t tapeJ1 = row[el];
                size_t tapeJ2 = col[el];
                hessLoop[tapeJ1][tapeJ2] = hessLoopFlat[el];
            }

            // store results in indexedLoopResults
            size_t hessElSize = info.indexedIndexedPositions.size() +
                    info.indexedTempPositions.size() +
                    info.nonIndexedIndexedPositions.size() +
                    info.nonIndexedNonIndexedPosition.size();
            map<size_t, map<pairss, std::vector<size_t> > >::const_iterator itJ1KJ2;

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

                indexedLoopResults[hessLE++] = createHessianContribution(handler, positions, hessLoop[tapeJ1].at(tapeJ2), *info.iterationIndexOp);
            }

            /**
             * indexed - non-indexed
             * -> done by  (non-indexed - indexed)
             */

            /**
             * indexed - temporary
             */
            map<pairss, set<size_t> >::const_iterator itEval;
            for (itEval = info.indexedTempEvals.begin(); itEval != info.indexedTempEvals.end(); ++itEval) {
                size_t tapeJ1 = itEval->first.first;
                size_t j2 = itEval->first.second;
                const set<size_t>& ks = itEval->second;

                const std::vector<HessianElement>& positions = info.indexedTempPositions[itEval->first];

                CGBase hessVal = Base(0);
                set<size_t>::const_iterator itz;
                for (itz = ks.begin(); itz != ks.end(); ++itz) {
                    size_t k = *itz;
                    size_t tapeK = lModel.getTempIndepIndexes(k)->tape;
                    hessVal += hessLoop[tapeJ1].at(tapeK) * dzDx[k][j2];
                }

                indexedLoopResults[hessLE++] = createHessianContribution(handler, positions, hessVal, *info.iterationIndexOp);
            }


            /*******************************************************************
             * non-indexed - indexed
             */
            for (it = info.nonIndexedIndexedPositions.begin(); it != info.nonIndexedIndexedPositions.end(); ++it) {
                size_t tapeJ1 = it->first.first;
                size_t tapeJ2 = it->first.second;
                const std::vector<HessianElement>& positions = it->second;

                indexedLoopResults[hessLE++] = createHessianContribution(handler, positions, hessLoop[tapeJ1].at(tapeJ2), *info.iterationIndexOp);
            }

            /*******************************************************************
             * temporary - indexed
             * 
             *      d f_i       .    d x_k1
             * d x_l2  d z_k1        d x_j1
             * 
             * -> done by  (indexed - temporary)
             */

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
                IndexPattern* pattern = new LinearIndexPattern(LoopModel<Base>::ITERATION_INDEX, 0, 0, 0, e);
                handler.manageLoopDependentIndexPattern(pattern);

                /**
                 * non-indexed - non-indexed
                 */
                CGBase hessVal = Base(0);

                if (info.nonIndexedNonIndexedEvals.find(orig) != info.nonIndexedNonIndexedEvals.end()) {
                    hessVal = hessLoop[posJ1->tape].at(posJ2->tape);
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
                        hessVal += hessLoop[posJ1->tape].at(tapeK) * dzDx[k][j2];
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
                        hessVal += hessLoop[tapeK].at(posJ2->tape) * dzDx[k1][j1];
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

                            tmp += hessLoop[tapeK1].at(tapeK2) * dzDx[k2][j2];
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
                    hessVal += dzDxx.at(j1).at(j2); // it is already the sum of ddz / dx_j1 dx_j2
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
            LoopEndOperationNode<Base>* loopEnd = createLoopEnd(handler, *info.loopStart, indexedLoopResults, indexesOps, lModel, assignOrAdd);

            std::vector<size_t> ninfo(1);
            std::vector<Argument<Base> > args(1);
            std::vector<size_t>::const_iterator itE;
            for (itE = lowerHessOrder.begin(); itE != lowerHessOrder.end(); ++itE) {
                // an additional alias variable is required so that each dependent variable can have its own ID
                size_t e = *itE;
                ninfo[0] = e;
                args[0] = Argument<Base>(*loopEnd);
                hess[e] = handler.createCG(new OperationNode<Base> (CGDependentRefOp, ninfo, args));
            }

            /**
             * move no-nindexed expressions outside loop
             */
            moveNonIndexedOutsideLoop(*info.loopStart, *loopEnd, LoopModel<Base>::ITERATION_INDEX);
        }

        return hess;
    }

    namespace loops {

        template<class Base>
        class HessianWithLoopsInfo {
        public:
            vector<std::set<size_t> > evalJacSparsity;
            vector<std::set<size_t> > evalHessSparsity;
            // (tapeJ1, tapeJ2) -> [positions]
            std::map<std::pair<size_t, size_t>, std::vector<HessianElement> > indexedIndexedPositions;
            // (tapeJ1, j2) -> [positions]
            std::map<std::pair<size_t, size_t>, std::vector<HessianElement> > indexedTempPositions;
            // (tapeJ1, j2) -> [k]
            std::map<std::pair<size_t, size_t>, std::set<size_t> > indexedTempEvals;
            // (tapeJ1(j1), tapeJ2) -> [positions]
            std::map<std::pair<size_t, size_t>, std::vector<HessianElement> > nonIndexedIndexedPositions;
            /**
             * (j1, j2) -> position
             */
            std::map<std::pair<size_t, size_t>, size_t> nonIndexedNonIndexedPosition;
            // [(j1, j2)]
            std::set<std::pair<size_t, size_t> > nonIndexedNonIndexedEvals;
            // (j2, j1) -> [k]
            std::map<std::pair<size_t, size_t>, std::set<size_t> > nonIndexedTempEvals;
            // (j1, j2) -> [k1]
            std::map<std::pair<size_t, size_t>, std::set<size_t> > tempNonIndexedEvals;
            // (j1, j2) -> k1 -> [k2]
            std::map<std::pair<size_t, size_t>, std::map<size_t, std::set<size_t> > > tempTempEvals;
            // (j1 ,j2) -> [k1]
            std::map<std::pair<size_t, size_t>, std::set<size_t> > nonLoopNonIndexedNonIndexed;

            LoopStartOperationNode<Base>* loopStart;
            IndexOperationNode<Base>* iterationIndexOp;
            vector<CG<Base> > x; // loop independent variables
            vector<CG<Base> > w;
            vector<std::map<size_t, CG<Base> > > dyiDzk;

            HessianWithLoopsInfo() :
                loopStart(NULL),
                iterationIndexOp(NULL) {

            }

            HessianWithLoopsInfo(CppAD::LoopModel<Base>& loop) :
                evalJacSparsity(loop.getTapeDependentCount()),
                evalHessSparsity(loop.getTapeIndependentCount()),
                loopStart(NULL),
                iterationIndexOp(NULL),
                dyiDzk(loop.getTapeDependentCount()) {

            }

        };

        template<class Base>
        std::pair<CG<Base>, IndexPattern*> createHessianContribution(CodeHandler<Base>& handler,
                                                                     const std::vector<HessianElement>& positions,
                                                                     const CG<Base>& ddfdxdx,
                                                                     IndexOperationNode<Base>& iterationIndexOp) {
            using namespace std;

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

                // generate the index pattern for the hessian compressed element
                IndexPattern* pattern = IndexPattern::detect(LoopModel<Base>::ITERATION_INDEX, locations.begin()->second);
                handler.manageLoopDependentIndexPattern(pattern);

                return make_pair(results.begin()->second, pattern);

            } else {
                /**
                 * must create a conditional element so that this 
                 * contribution to the hessian is only evaluated at the
                 * relevant iterations
                 */
                const std::vector<size_t> ninfo;
                OperationNode<Base>* ifBranch;
                std::vector<Argument<Base> > nextBranchArgs;
                set<size_t> usedIter;

                map<size_t, map<size_t, size_t> >::const_iterator countIt;
                for (countIt = locations.begin(); countIt != locations.end(); ++countIt) {
                    size_t count = countIt->first;
                    const map<size_t, size_t>& locationsC = countIt->second;
                    size_t iterCount = locationsC.size();

                    if (usedIter.size() + iterCount == positions.size()) {
                        // all other iterations
                        ifBranch = new OperationNode<Base>(CGElseOp, ninfo, nextBranchArgs);
                        handler.manageOperationNodeMemory(ifBranch);
                    } else {
                        std::set<size_t> iterations;
                        mapKeys(locationsC, iterations);

                        // depends on the iterations indexes
                        OperationNode<Base>* cond = createIndexConditionExpression<Base>(iterations, usedIter, positions.size() - 1, iterationIndexOp);
                        handler.manageOperationNodeMemory(cond);
                        //new OperationNode<Base>(CGCondExprOp, ninfo, args);
                        nextBranchArgs.insert(nextBranchArgs.begin(), Argument<Base>(*cond));
                        ifBranch = new OperationNode<Base>(countIt == locations.begin() ? CGStartIfOp : CGElseIfOp, ninfo, nextBranchArgs);
                        handler.manageOperationNodeMemory(ifBranch);

                        usedIter.insert(iterations.begin(), iterations.end());
                    }

                    nextBranchArgs.clear();

                    IndexPattern* pattern = IndexPattern::detect(LoopModel<Base>::ITERATION_INDEX, locationsC);
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
                    nextBranchArgs.push_back(Argument<Base>(*ifAssign));
                }

                OperationNode<Base>* endif = new OperationNode<Base>(CGEndIfOp, ninfo, nextBranchArgs);

                IndexPattern* pattern = NULL;
                return make_pair(handler.createCG(endif), pattern);
            }

        }

        template<class Base>
        OperationNode<Base>* createIndexConditionExpression(const std::set<size_t>& iterations,
                                                            const std::set<size_t>& usedIter,
                                                            size_t maxIter,
                                                            IndexOperationNode<Base>& iterationIndexOp) {
            assert(!iterations.empty());

            std::map<size_t, bool> allIters;
            for (std::set<size_t>::const_iterator it = usedIter.begin(); it != usedIter.end(); ++it) {
                allIters[*it] = false;
            }
            for (std::set<size_t>::const_iterator it = iterations.begin(); it != iterations.end(); ++it) {
                allIters[*it] = true;
            }

            std::vector<size_t> info;
            info.reserve(iterations.size() / 2 + 2);

            std::map<size_t, bool>::const_iterator it = allIters.begin();
            while (it != allIters.end()) {
                std::map<size_t, bool> ::const_iterator min = it;
                std::map<size_t, bool>::const_iterator max = it;
                std::map<size_t, bool>::const_iterator minNew = allIters.end();
                std::map<size_t, bool>::const_iterator maxNew = allIters.end();
                if (it->second) {
                    minNew = it;
                    maxNew = it;
                }

                for (++it; it != allIters.end(); ++it) {
                    if (it->first != max->first + 1) {
                        break;
                    }

                    max = it;
                    if (it->second) {
                        if (minNew == allIters.end())
                            minNew = it;
                        maxNew = it;
                    }
                }

                if (minNew != allIters.end()) {
                    // contains elements from the current iteration set
                    if (maxNew->first == minNew->first) {
                        // only one element
                        info.push_back(minNew->first);
                        info.push_back(maxNew->first);
                    } else {
                        //several elements
                        if (min->first == 0)
                            info.push_back(min->first);
                        else
                            info.push_back(minNew->first);

                        if (max->first == maxIter)
                            info.push_back(std::numeric_limits<size_t>::max());
                        else
                            info.push_back(maxNew->first);
                    }
                }
            }

            std::vector<Argument<Base> > args(1);
            args[0] = Argument<Base>(iterationIndexOp);

            return new OperationNode<Base>(CGIndexCondExprOp, info, args);
        }
    }
}

#endif
