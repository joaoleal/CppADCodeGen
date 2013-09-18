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

        inline void evalHessianSparsity() {
            if (!hessSparsity_) {
                hessTapeSparsity_ = hessianSparsitySet<vector<std::set<size_t> >, CGB>(*fun_);
                hessSparsity_ = true;
            }
        }

        inline const vector<std::set<size_t> >& getHessianSparsity() const {
            return hessTapeSparsity_;
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