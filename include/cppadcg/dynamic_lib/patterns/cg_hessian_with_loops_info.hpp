#ifndef CPPAD_CG_HESSIAN_WITH_LOOPS_INFO_INCLUDED
#define CPPAD_CG_HESSIAN_WITH_LOOPS_INFO_INCLUDED
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

        template<class Base>
        class HessianWithLoopsInfo {
        public:
            CppAD::LoopModel<Base>* model;
            vector<std::set<size_t> > evalJacSparsity;
            vector<std::set<size_t> > evalHessSparsity;
            // (tapeJ1, tapeJ2) -> [positions]
            std::map<pairss, std::vector<HessianElement> > indexedIndexedPositions;
            // (tapeJ1, tapeJ2(j2)) -> [positions]
            std::map<pairss, std::vector<HessianElement> > indexedNonIndexedPositions;
            // (tapeJ1, j2) -> [positions]
            std::map<pairss, std::vector<HessianElement> > indexedTempPositions;

            // (tapeJ1(j1), tapeJ2) -> [positions]
            std::map<pairss, std::vector<HessianElement> > nonIndexedIndexedPositions;
            //(j1, j2) -> position
            std::map<pairss, size_t> nonIndexedNonIndexedPosition;

            // (j1, tapeJ2) -> [positions]
            std::map<pairss, std::vector<HessianElement> > tempIndexedPositions;

            // (tapeJ1, j2) -> [k]
            std::map<pairss, std::set<size_t> > indexedTempEvals;
            // [(j1, j2)]
            std::set<pairss> nonIndexedNonIndexedEvals;
            // (j2, j1) -> [k]
            std::map<pairss, std::set<size_t> > nonIndexedTempEvals;
            // (j1, j2) -> [k1]
            std::map<pairss, std::set<size_t> > tempNonIndexedEvals;
            // (j1, j2) -> k1 -> [k2]
            std::map<pairss, std::map<size_t, std::set<size_t> > > tempTempEvals;
            // (j1 ,j2) -> [k1]
            std::map<pairss, std::set<size_t> > nonLoopNonIndexedNonIndexed;

            LoopStartOperationNode<Base>* loopStart;
            LoopEndOperationNode<Base>* loopEnd;
            IndexOperationNode<Base>* iterationIndexOp;
            vector<CG<Base> > x; // loop independent variables
            vector<CG<Base> > w;
            vector<std::map<size_t, CG<Base> > > dyiDzk;

            vector<std::map<size_t, CG<Base> > > hess;

            vector<std::set<size_t> > noLoopEvalHessTempsSparsity;
            std::map<size_t, std::map<size_t, CG<Base> > > dzDxx;

            // if-else branches
            vector<IfElseInfo<Base> > ifElses;

            HessianWithLoopsInfo() :
                model(NULL),
                loopStart(NULL),
                loopEnd(NULL),
                iterationIndexOp(NULL) {

            }

            HessianWithLoopsInfo(CppAD::LoopModel<Base>& loop) :
                model(&loop),
                evalJacSparsity(loop.getTapeDependentCount()),
                evalHessSparsity(loop.getTapeIndependentCount()),
                loopStart(NULL),
                loopEnd(NULL),
                iterationIndexOp(NULL) {

            }

            inline void evalLoopModelJacobian() {
                ADFun<CG<Base> >& fun = model->getTape();
                const vector<std::set<size_t> >& jacTapeSparsity = model->getJacobianSparsity();

                //printSparsityPattern(jacTapeSparsity, "jac - loop");
                //printSparsityPattern(info.evalJacSparsity, "jac - loop -eval");

                /**
                 * evaluate loop model jacobian
                 */
                std::vector<size_t> row, col;
                extra::generateSparsityIndexes(evalJacSparsity, row, col);
                if (row.size() > 0) {
                    vector<CG<Base> > jacLoop(row.size());

                    CppAD::sparse_jacobian_work work; // temporary structure for CppAD
                    if (extra::estimateBestJacobianADMode(row, col)) {
                        fun.SparseJacobianForward(x, jacTapeSparsity, row, col, jacLoop, work);
                    } else {
                        fun.SparseJacobianReverse(x, jacTapeSparsity, row, col, jacLoop, work);
                    }

                    // save/organize results
                    dyiDzk.resize(model->getTapeDependentCount());
                    for (size_t el = 0; el < jacLoop.size(); el++) {
                        size_t tapeI = row[el];
                        size_t tapeJ = col[el];
                        dyiDzk[tapeI][tapeJ] = jacLoop[el];
                    }
                }

            }

            inline void evalLoopModelJacobianHessian() {
                using namespace CppAD::extra;
                
                ADFun<CG<Base> >& fun = model->getTape();

                std::vector<size_t> jacRow, jacCol;
                generateSparsityIndexes(evalJacSparsity, jacRow, jacCol);
                vector<CG<Base> > jacLoop(jacRow.size());

                std::vector<size_t> hesRow, hesCol;
                generateSparsityIndexes(evalHessSparsity, hesRow, hesCol);
                vector<CG<Base> > hessLoopFlat(hesRow.size());

                if (!jacRow.empty() || !hesRow.empty()) {
                    vector<CG<Base> > y(model->getTapeDependentCount());

                    SparseForjacHessianWork work;
                    sparseForJacHessian(fun, x, w,
                                        y,
                                        model->getJacobianSparsity(),
                                        jacRow, jacCol, jacLoop,
                                        model->getHessianSparsity(),
                                        hesRow, hesCol, hessLoopFlat,
                                        work);

                    // save/organize results
                    // Jacobian
                    dyiDzk.resize(model->getTapeDependentCount());
                    for (size_t el = 0; el < jacLoop.size(); el++) {
                        size_t tapeI = jacRow[el];
                        size_t tapeJ = jacCol[el];
                        dyiDzk[tapeI][tapeJ] = jacLoop[el];
                    }

                    // save non-indexed hessian elements
                    // Hessian
                    hess.resize(fun.Domain());
                    for (size_t el = 0; el < hesRow.size(); el++) {
                        size_t tapeJ1 = hesRow[el];
                        size_t tapeJ2 = hesCol[el];
                        hess[tapeJ1][tapeJ2] = hessLoopFlat[el];
                    }
                }
            }

        };

    }
}

#endif