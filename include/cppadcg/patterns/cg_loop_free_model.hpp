#ifndef CPPAD_CG_LOOP_FREE_MODEL_INCLUDED
#define CPPAD_CG_LOOP_FREE_MODEL_INCLUDED
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

    /**
     * Altered model without the loop equations and with extra dependents
     * for the non-indexed temporary variables used by loops
     * 
     * @author Joao Leal
     */
    template <class Base>
    class LoopFreeModel {
    public:
        typedef CppAD::CG<Base> CGB;
        typedef Argument<Base> Arg;
    protected:
        /**
         * The tape
         */
        ADFun<CGB> * const fun_;
        /**
         * The dependent variables in this tape to their original indexes
         */
        std::vector<size_t> dependentIndexes_;
        std::map<size_t, size_t> dependentOrig2Local;
        /**
         * Jacobian sparsity pattern of the tape
         */
        vector<std::set<size_t> > jacTapeSparsity_;
        bool jacSparsity_;
        /**
         * Hessian sparsity pattern for equations used to determine the 
         * temporaries (ignores the the original model equations)
         */
        vector<std::set<size_t> > hessTapeTempSparsity_;
        /**
         * Hessian sparsity pattern for the original model equations in the tape
         * (ignores the equations for the temporaries)
         */
        vector<std::set<size_t> > hessTapeOrigEqSparsity_;
        // whether or not the hessian sparsities have been evaluated
        bool hessSparsity_;
    public:

        /**
         * Creates a model for the non-indexed operations
         * 
         * @param fun
         * @param dependentOrigIndexes
         */
        LoopFreeModel(ADFun<CGB>* fun,
                      const std::vector<size_t>& dependentOrigIndexes) :
            fun_(fun),
            dependentIndexes_(dependentOrigIndexes),
            jacSparsity_(false),
            hessSparsity_(false) {
            CPPADCG_ASSERT_KNOWN(fun != NULL, "fun cannot be NULL");
            CPPADCG_ASSERT_KNOWN(dependentOrigIndexes.size() <= fun->Range(), "invalid size");

            for (size_t il = 0; il < dependentIndexes_.size(); il++)
                dependentOrig2Local[dependentIndexes_[il]] = il;
        }

        inline ADFun<CGB>& getTape() const {
            return *fun_;
        }

        inline size_t getTapeDependentCount() const {
            return fun_->Range();
        }

        inline size_t getTemporaryDependentCount() const {
            return fun_->Range() - dependentIndexes_.size();
        }

        inline size_t getTapeIndependentCount() const {
            return fun_->Domain();
        }

        /**
         * Provides the dependent variables indexes present in the original
         * model
         */
        inline const std::vector<size_t>& getOrigDependentIndexes() const {
            return dependentIndexes_;
        }

        inline size_t getLocalDependentIndex(size_t origI) const {
            return dependentOrig2Local.at(origI);
        }

        inline void evalJacobianSparsity() {
            if (!jacSparsity_) {
                jacTapeSparsity_ = jacobianSparsitySet<vector<std::set<size_t> >, CGB>(*fun_);
                jacSparsity_ = true;
            }
        }

        inline const vector<std::set<size_t> >& getJacobianSparsity() const {
            return jacTapeSparsity_;
        }

        inline void evalHessianSparsity() {
            if (!hessSparsity_) {
                size_t mo = dependentIndexes_.size();
                size_t m = fun_->Range();
                size_t n = fun_->Domain();

                // hessian for the original equations
                std::set<size_t> eqs;
                if (mo != 0) {
                    for (size_t i = 0; i < mo; i++)
                        eqs.insert(eqs.end(), i);

                    hessTapeOrigEqSparsity_ = hessianSparsitySet<vector<std::set<size_t> >, CGB>(*fun_, eqs);
                }

                // hessian for the temporary variable equations
                if (m != mo) {
                    eqs.clear();
                    for (size_t i = mo; i < m; i++)
                        eqs.insert(eqs.end(), i);
                    hessTapeTempSparsity_ = hessianSparsitySet<vector<std::set<size_t> >, CGB>(*fun_, eqs);
                } else {
                    hessTapeTempSparsity_.resize(n);
                }

                hessSparsity_ = true;
            }
        }

        inline const vector<std::set<size_t> >& getHessianTempEqsSparsity() const {
            assert(hessSparsity_);
            return hessTapeTempSparsity_;
        }

        inline const vector<std::set<size_t> >& getHessianOrigEqsSparsity() const {
            assert(hessSparsity_);
            return hessTapeOrigEqSparsity_;
        }

        /**
         * Determines the hessian for the temporary variables only used by
         * each loop
         * 
         * @param loopHessInfo loops
         * @param x the independent variables
         */
        inline std::map<size_t, std::map<size_t, CGB> > calculateJacobianHessianUsedByLoops(std::map<LoopModel<Base>*, loops::HessianWithLoopsInfo<Base> >& loopHessInfo,
                                                                                            const vector<CGB>& x,
                                                                                            vector<CGB>& temps,
                                                                                            const vector<std::set<size_t> >& noLoopEvalJacSparsity) {
            using namespace std;
            using namespace CppAD::loops;

            assert(hessSparsity_); // check that the sparsities have been evaluated

            size_t mo = dependentIndexes_.size();
            size_t m = getTapeDependentCount();
            size_t n = fun_->Domain();

            vector<vector<CGB> > vwNoLoop(loopHessInfo.size());
            vector<vector<CGB> > vhessNoLoop(loopHessInfo.size());

            // jacobian for temporaries
            std::vector<size_t> jacRow, jacCol;
            generateSparsityIndexes(noLoopEvalJacSparsity, jacRow, jacCol);

            // jacobian for equations outside loops
            vector<CGB> jacNoLoop(jacRow.size());

            /**
             * hessian - temporary variables
             */
            vector<std::set<size_t> > noLoopEvalHessTempsSparsity(n);

            typename map<LoopModel<Base>*, HessianWithLoopsInfo<Base> >::iterator itLoop2Info;
            for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info) {
                HessianWithLoopsInfo<Base>& info = itLoop2Info->second;

                addMatrixSparsity(info.noLoopEvalHessTempsSparsity, noLoopEvalHessTempsSparsity);
            }
            std::vector<size_t> hesRow, hesCol;
            generateSparsityIndexes(noLoopEvalHessTempsSparsity, hesRow, hesCol);

            size_t l = 0;
            for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info, l++) {
                LoopModel<Base>* loop = itLoop2Info->first;
                HessianWithLoopsInfo<Base>& info = itLoop2Info->second;

                vector<CGB>& hessNoLoop = vhessNoLoop[l];
                hessNoLoop.resize(hesRow.size());

                vector<CGB>& wNoLoop = vwNoLoop[l];
                wNoLoop.resize(m);
                for (size_t inl = 0; inl < mo; inl++) {
                    wNoLoop[inl] = Base(0);
                }

                for (size_t inl = mo; inl < m; inl++) {
                    size_t k = inl - mo;
                    const LoopPosition* posK = loop->getTempIndepIndexes(k);

                    if (posK != NULL) {
                        for (size_t i = 0; i < info.dyiDzk.size(); i++) {
                            const map<size_t, CGB>& row = info.dyiDzk[i];
                            typename map<size_t, CGB>::const_iterator itCol = row.find(posK->tape);
                            if (itCol != row.end()) {
                                const CGB& val = itCol->second;
                                wNoLoop[inl] += val * info.w[i];
                            }
                        }
                    }
                }
            }

            SparseForjacHessianWork work;
            sparseForJacHessian(*fun_, x, vwNoLoop,
                                temps,
                                getJacobianSparsity(),
                                jacRow, jacCol, jacNoLoop,
                                hessTapeTempSparsity_,
                                hesRow, hesCol, vhessNoLoop,
                                work);

            // save Jacobian
            std::map<size_t, std::map<size_t, CGB> > dzDx;
            for (size_t el = 0; el < jacRow.size(); el++) {
                size_t inl = jacRow[el];
                size_t j = jacCol[el];
                assert(inl >= mo);

                // dz_k/dx_v (for temporary variable)
                size_t k = inl - mo;
                dzDx[k][j] = jacNoLoop[el];
            }

            // save Hessian
            l = 0;
            for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info, l++) {
                HessianWithLoopsInfo<Base>& info = itLoop2Info->second;
                vector<CGB>& hessNoLoop = vhessNoLoop[l];

                for (size_t el = 0; el < hesRow.size(); el++) {
                    size_t j1 = hesRow[el];
                    size_t j2 = hesCol[el];
                    info.dzDxx[j1][j2] = hessNoLoop[el];
                }
            }

            return dzDx;
        }

        inline void calculateHessian4OrignalEquations(const vector<CGB>& x,
                                                      const vector<CGB>& w,
                                                      const vector<std::set<size_t> >& noLoopEvalHessSparsity,
                                                      const vector<std::map<size_t, std::set<size_t> > >& noLoopEvalHessLocations,
                                                      vector<CGB>& hess) {
            using namespace std;
            using namespace CppAD::loops;

            assert(hessSparsity_); // check that the sparsities have been evaluated

            vector<CGB> wNoLoop(getTapeDependentCount());
            vector<CGB> hessNoLoop;

            /**
             * hessian - original equations
             */
            std::vector<size_t> row, col;
            generateSparsityIndexes(noLoopEvalHessSparsity, row, col);

            if (row.size() > 0) {
                hessNoLoop.resize(row.size());

                for (size_t inl = 0; inl < dependentIndexes_.size(); inl++) {
                    wNoLoop[inl] = w[dependentIndexes_[inl]];
                }

                CppAD::sparse_hessian_work work; // temporary structure for CPPAD
                fun_->SparseHessian(x, wNoLoop, hessTapeOrigEqSparsity_, row, col, hessNoLoop, work);

                // save non-indexed hessian elements
                for (size_t el = 0; el < row.size(); el++) {
                    size_t j1 = row[el];
                    size_t j2 = col[el];
                    const set<size_t>& locations = noLoopEvalHessLocations[j1].at(j2);
                    for (set<size_t>::const_iterator itE = locations.begin(); itE != locations.end(); ++itE)
                        hess[*itE] = hessNoLoop[el];
                }
            }
        }

        virtual ~LoopFreeModel() {
            delete fun_;
        }

    private:
        LoopFreeModel(const LoopFreeModel<Base>&); // not implemented

        LoopFreeModel& operator=(const LoopFreeModel<Base>&); // not implemented
    };

}

#endif