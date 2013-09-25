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
         * Hessian sparsity pattern of the tape
         */
        vector<std::set<size_t> > hessTapeSparsity_;
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

        inline std::map<size_t, std::map<size_t, CGB> > evaluateJacobian4Temporaries(const vector<CGB>& x,
                                                                                     const vector<std::set<size_t> >& noLoopEvalJacSparsity) {
            using namespace std;

            size_t nonIndexdedEqSize = dependentIndexes_.size();

            // jacobian for equations outside loops
            vector<CGB> jacNoLoop;
            // jacobian for temporaries
            map<size_t, map<size_t, CGB> > dzDx;

            std::vector<size_t> row, col;
            generateSparsityIndexes(noLoopEvalJacSparsity, row, col);

            if (row.size() > 0) {
                jacNoLoop.resize(row.size());

                //printSparsityPattern(_funNoLoops->getJacobianSparsity(), "jacobian No Loops");
                //printSparsityPattern(noLoopEvalJacSparsity, "jacobian No Loops - eval");

                CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
                if (estimateBestJacobianADMode(row, col)) {
                    fun_->SparseJacobianForward(x, getJacobianSparsity(), row, col, jacNoLoop, work);
                } else {
                    fun_->SparseJacobianReverse(x, getJacobianSparsity(), row, col, jacNoLoop, work);
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

            return dzDx;
        }

        inline void evalHessianSparsity() {
            if (!hessSparsity_) {
                hessTapeSparsity_ = hessianSparsitySet<vector<std::set<size_t> >, CGB>(*fun_);
                hessSparsity_ = true;
            }
        }

        inline const vector<std::set<size_t> >& getHessianSparsity() const {
            return hessTapeSparsity_;
        }

        /**
         * Determines the hessian for the temporary variables only used by
         * each loop
         * 
         * @param loopHessInfo loops
         * @param x the independent variables
         */
        inline void calculateHessianUsedByLoops(std::map<LoopModel<Base>*, loops::HessianWithLoopsInfo<Base> >& loopHessInfo,
                                                const vector<CGB>& x) {
            using namespace std;
            using namespace CppAD::loops;

            vector<CGB> wNoLoop(getTapeDependentCount());

            vector<CGB> hessNoLoop;

            /**
             * hessian - temporary variables
             */
            std::vector<size_t> row, col;

            typename map<LoopModel<Base>*, HessianWithLoopsInfo<Base> >::iterator itLoop2Info;
            for (itLoop2Info = loopHessInfo.begin(); itLoop2Info != loopHessInfo.end(); ++itLoop2Info) {
                LoopModel<Base>* loop = itLoop2Info->first;
                HessianWithLoopsInfo<Base>& info = itLoop2Info->second;

                generateSparsityIndexes(info.noLoopEvalHessTempsSparsity, row, col);

                if (row.empty())
                    continue;

                hessNoLoop.resize(row.size());

                for (size_t inl = 0; inl < dependentIndexes_.size(); inl++) {
                    wNoLoop[inl] = Base(0);
                }

                for (size_t inl = dependentIndexes_.size(); inl < wNoLoop.size(); inl++) {
                    size_t k = inl - dependentIndexes_.size();
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

                CppAD::sparse_hessian_work workTemps;
                fun_->SparseHessian(x, wNoLoop, getHessianSparsity(), row, col, hessNoLoop, workTemps);

                // save hessian
                for (size_t el = 0; el < row.size(); el++) {
                    size_t j1 = row[el];
                    size_t j2 = col[el];
                    info.dzDxx[j1][j2] = hessNoLoop[el];
                }
            }
        }

        inline void calculateHessian4OrignalEquations(const vector<CGB>& x,
                                                      const vector<CGB>& w,
                                                      const vector<std::set<size_t> >& noLoopEvalHessSparsity,
                                                      const vector<std::map<size_t, std::set<size_t> > >& noLoopEvalHessLocations,
                                                      vector<CGB>& hess) {
            using namespace std;
            using namespace CppAD::loops;

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
                fun_->SparseHessian(x, wNoLoop, getHessianSparsity(), row, col, hessNoLoop, work);

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