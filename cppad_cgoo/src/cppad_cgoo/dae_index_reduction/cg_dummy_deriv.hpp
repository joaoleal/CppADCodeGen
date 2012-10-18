#ifndef CPPAD_CG_DUMMY_DERIV_INCLUDED
#define	CPPAD_CG_DUMMY_DERIV_INCLUDED

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/NonLinearOptimization>

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
namespace CppAD {

    template<class Base>
    bool sortVnodesByOrder(Vnode<Base>* i, Vnode<Base>* j) {
        return (i->order() > j->order());
    }

    /**
     * Dummy derivatives DAE index reduction algorithm
     */
    template<class Base>
    class DummyDerivatives : public Plantelides<Base> {
        typedef Eigen::Matrix<Base, Eigen::Dynamic, 1 > VectorB;
        typedef Eigen::Matrix<std::complex<Base>, Eigen::Dynamic, 1 > VectorCB;
        typedef Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic> MatrixB;
    protected:
        // typical values;
        std::vector<Base> x_;
        // normalization constants
        std::vector<Base> norm_;
        // jacobian sparsity pattern of the reduced system
        std::vector<bool> jacSparsity_;
        // the initial index of time derivatives
        size_t diffVarStart_;
        // the initial index of the differentiated equations
        size_t diffEqStart_;
        /** normalized Jacobian of the index one system's  differentiated
         *  equations relative to the time derivatives
         */
        Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic> jacobian_;
        /**
         * Dummy derivatives
         */
        std::set<Vnode<Base>* > dummyD_;
    public:

        DummyDerivatives(ADFun<CG<Base> >* fun,
                         const std::vector<bool>& eqDifferentialInfo,
                         const std::vector<bool>& varInfo,
                         const std::vector<Base>& x,
                         const std::vector<Base>& norm) :
            Plantelides<Base>(fun, eqDifferentialInfo, varInfo),
            x_(x),
            norm_(norm),
            diffVarStart_(0),
            diffEqStart_(this->eqDifferentialInfo_.size()) {

            typename std::vector<Vnode<Base>*> ::const_iterator j;
            for (j = this->vnodes_.begin(); j != this->vnodes_.end(); ++j) {
                Vnode<Base>* jj = *j;
                if (jj->derivativeOf() != NULL) {
                    diffVarStart_ = jj->index();
                    break;
                }
            }
        }

        virtual inline void reduceIndex() {
            Plantelides<Base>::reduceIndex();

            assert(this->reducedFun_ != NULL);

            //solveDAESystem();

            determineJacobian();


            // variables of interest
            std::vector<Vnode<Base>* > vars;
            vars.reserve(this->vnodes_.size() - diffVarStart_);
            typename std::vector<Vnode<Base>* >::const_reverse_iterator rj;
            for (rj = this->vnodes_.rbegin(); rj != this->vnodes_.rend(); ++rj) {
                Vnode<Base>* jj = *rj;
                if (jj->derivativeOf() != NULL && jj->derivative() == NULL) {
                    vars.push_back(jj);
                }
            }

            // should be already fairly sorted, but sort anyway
            std::sort(vars.begin(), vars.end(), sortVnodesByOrder<Base>);

            // equations of interest
            typename std::vector<Enode<Base>* >::const_reverse_iterator ri;
            std::vector<Enode<Base>* > eqs;
            eqs.reserve(this->enodes_.size() - diffEqStart_);
            for (ri = this->enodes_.rbegin(); ri != this->enodes_.rend(); ++ri) {
                Enode<Base>* ii = *ri;
                if (ii->derivativeOf() != NULL && ii->derivative() == NULL) {
                    eqs.push_back(ii);
                }
            }


            MatrixB workJac;

            while (true) {

#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << "# equation selection: ";
                for (size_t i = 0; i < eqs.size(); i++)
                    std::cout << *eqs[i] << "; ";
                std::cout << "\n";

                std::cout << "# variable selection: ";
                for (size_t j = 0; j < vars.size(); j++)
                    std::cout << *vars[j] << "; ";
                std::cout << "\n";
#endif

                // create the jacobian for the selected variables and equations
                workJac.resize(eqs.size(), vars.size());
                for (size_t i = 0; i < eqs.size(); i++) {
                    Enode<Base>* ii = eqs[i];
                    for (size_t j = 0; j < vars.size(); j++) {
                        Vnode<Base>* jj = vars[j];
                        workJac(i, j) = jacobian_(ii->index() - diffEqStart_, jj->index() - diffVarStart_);
                    }
                }

                // Exploit the current equations for elimination of candidates
                selectDummyDerivatives(eqs, vars, workJac);

                /**
                 * Consider all of the current equations that are
                 * differentiated versions of the original ones.
                 * Collect their predecessors and let them be the
                 * current equations.
                 */
                std::vector<Enode<Base>* > newEqs;
                newEqs.reserve(eqs.size());

                typename std::vector<Enode<Base>* >::const_iterator i;
                for (i = eqs.begin(); i != eqs.end(); ++i) {
                    Enode<Base>* ii = (*i)->derivativeOf();
                    if (ii != NULL && ii->derivativeOf() != NULL) {
                        newEqs.push_back(ii);
                    }
                }
                eqs.swap(newEqs);

                if (eqs.empty()) {
                    break;
                }

                /**
                 * Consider all current unknowns that are at least of
                 * order one. Collect their predecessors of one order
                 * less and let them be the current candidates for
                 * elimination.
                 */
                std::vector<Vnode<Base>* > varsNew;
                varsNew.reserve(vars.size());
                typename std::vector<Vnode<Base>* >::const_iterator j;
                for (j = vars.begin(); j != vars.end(); ++j) {
                    Vnode<Base>* v = (*j)->derivativeOf();
                    if (v != NULL && v->derivativeOf() != NULL) {
                        varsNew.push_back(v);
                    }
                }
                vars.swap(varsNew);

                size_t newDummyCount = vars.size() - eqs.size();
            }
        }

    protected:

        inline void solveDAESystem() {
            Functor dae(this);
            Eigen::LevenbergMarquardt<Functor> lm(dae);

            size_t size = dae.inputs(); // number of equations and variables

            VectorB x(size);
            for (size_t j = 0, pos = 0; j< this->eqDifferentialInfo_.size(); j++) {
                if (this->eqDifferentialInfo_[j]) {
                    x(pos++) = x_[j];
                }
            }

            int info = lm.minimize(x);
        }

        /**
         * Determines the Jacobian relative to the differential variables
         * (e.g. dxdt)
         */
        inline void determineJacobian() {
            const size_t n = this->reducedFun_->Domain();
            const size_t m = this->reducedFun_->Range();

            jacSparsity_ = jacobianReverseSparsity(*this->reducedFun_);

            std::vector<size_t> row, col;
            row.reserve((this->vnodes_.size() - diffVarStart_) * (m - diffEqStart_));
            col.reserve(row.capacity());

            for (size_t i = diffEqStart_; i < m; i++) {
                for (size_t j = diffVarStart_; j < n; j++) {
                    if (jacSparsity_[i * n + j]) {
                        row.push_back(i);
                        col.push_back(j);
                    }
                }
            }

            std::vector<CG<Base> > jac(row.size());

            std::vector<CG<Base> > indep(n);
            std::copy(x_.begin(), x_.end(), indep.begin());
            std::fill(indep.begin() + x_.size(), indep.end(), 0);

            CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
            this->reducedFun_->SparseJacobianReverse(indep, jacSparsity_,
                                                     row, col, jac, work);

            // resize and zero matrix
            jacobian_.setZero(m - diffEqStart_, n - diffVarStart_);

            // right hand side of the equations
            for (size_t e = 0; e < jac.size(); e++) {
                Vnode<Base>* v = this->vnodes_[col[e]]->originalVariable();
                Enode<Base>* eq = this->enodes_[row[e]]->originalEquation();

                Base normVal = jac[e].getParameterValue() * norm_[v->index()]; // normalized jacobian value
                if (this->eqDifferentialInfo_[eq->index()]) {
                    normVal /= norm_[eq->index()];
                }

                jacobian_(row[e] - diffEqStart_, col[e] - diffVarStart_) = normVal;
            }

            // left hand side of the equations
            for (size_t i = 0; i < diffEqStart_; i++) {
                if (this->eqDifferentialInfo_[i]) {
                    Vnode<Base>* vDiff = this->vnodes_[i]->derivative();
                    Enode<Base>* eq = this->enodes_[i];

                    while (eq->derivative() != NULL) {
                        assert(vDiff != NULL);

                        eq = eq->derivative();
                        vDiff = vDiff->derivative();

                        size_t r = eq->index() - diffEqStart_;
                        size_t c = vDiff->index() - diffVarStart_;
                        jacobian_(r, c) = -1.0 / norm_[i];
                    }
                }
            }

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "partial jacobian:\n" << jacobian_ << "\n\n";
#endif
        }

        inline void selectDummyDerivatives(const std::vector<Enode<Base>* >& eqs,
                                           const std::vector<Vnode<Base>* >& vars,
                                           MatrixB& subsetJac) {

            if (eqs.size() == vars.size()) {
                dummyD_.insert(vars.begin(), vars.end());
#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << "# new dummy derivatives: ";
                for (size_t j = 0; j < vars.size(); j++)
                    std::cout << *vars[j] << "; ";
                std::cout << " \n";
#endif
                return;
            }

            subsetJac.resize(eqs.size(), vars.size());
            for (size_t i = 0; i < eqs.size(); i++) {
                Enode<Base>* ii = eqs[i];
                for (size_t j = 0; j < vars.size(); j++) {
                    Vnode<Base>* jj = vars[j];
                    subsetJac(i, j) = jacobian_(ii->index() - diffEqStart_, jj->index() - diffVarStart_);
                }
            }

            // number of columns/variables to remove (the remaining will be dummy derivatives)
            std::vector<size_t> cols2keep(eqs.size());
            for (size_t c = 0; c < eqs.size(); c++) {
                cols2keep[c] = c;
            }

            MatrixB workJac(eqs.size(), eqs.size());

            /**
             * Brute force approach!!!
             */
            std::vector<size_t> bestCols2keep = cols2keep;
            Base bestCond = std::numeric_limits<Base>::max();
            size_t bestTotalOrder = 0;

            while (true) {

#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << " ## column selection: ";
                for (size_t s = 0; s < cols2keep.size(); s++)
                    std::cout << cols2keep[s] << " ";
                std::cout << " \n";
#endif

                for (size_t i = 0; i < eqs.size(); i++) {
                    for (size_t j = 0; j < cols2keep.size(); j++) {
                        workJac(i, j) = subsetJac(i, cols2keep[j]);
                    }
                }
#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << "    current jac:\n" << workJac << "\n";
#endif

                Base cond = evalCondition(workJac);

#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << "    condition: " << cond << "\n";
#endif

                if (cond == cond) {
                    // not NaN
                    size_t totalOrd = 0;
                    for (size_t j = 0; j < cols2keep.size(); j++) {
                        totalOrd += vars[cols2keep[j]]->order();
                    }
                    if (cond <= bestCond) {
                        if (totalOrd >= bestTotalOrder) {
                            bestTotalOrder = totalOrd;
                            bestCond = cond;
                            bestCols2keep = cols2keep;
                        }
                    }
                }

                /**
                 * determine the next set of columns
                 */
                if (cols2keep.back() == vars.size() - 1) {
                    if (cols2keep[0] == vars.size() - cols2keep.size())
                        break; // end of combinations

                    for (size_t cc = 1; cc < cols2keep.size(); cc++) {
                        if (cols2keep[cc] == vars.size() - (cols2keep.size() - cc)) {
                            cols2keep[cc - 1]++;
                            for (size_t cc2 = cc; cc2 < cols2keep.size(); cc2++) {
                                cols2keep[cc2] = cols2keep[cc2 - 1] + 1;
                            }
                            break;
                        }
                    }
                } else {
                    cols2keep.back()++;
                }
            };

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "# new dummy derivatives (condition = " << bestCond << "): ";
            for (size_t c = 0; c < bestCols2keep.size(); c++)
                std::cout << *vars[bestCols2keep[c]] << "; ";
            std::cout << " \n";
#endif

            for (size_t c = 0; c < bestCols2keep.size(); c++) {
                dummyD_.insert(vars[bestCols2keep[c]]);
            }

        }

        /**
         * 
         */
        struct Functor {
            const DummyDerivatives<Base> * const dummyDer_;
            ADFun<Base>* reducedFunB_;
            std::vector<Base> normdep_;
            std::vector<Base> normindep_;
            std::vector<Base> jac_; // Jacobian
            std::vector<size_t> row_; // Jacobian row indexes
            std::vector<size_t> col_; // Jacobian column indexes
            std::vector< std::set<size_t> > jac_sparsity_; // Jacobian column indexes
            CppAD::sparse_jacobian_work work_; // temporary structure for CPPAD

            Functor(DummyDerivatives<Base>* dummyDer) :
                dummyDer_(dummyDer),
                normdep_(dummyDer_->reducedFun_->Range(), 1.0),
                normindep_(dummyDer_->reducedFun_->Range(), 1.0) {

                /**
                 * get rid of the CG encapsulation
                 */
                CodeHandler<Base> handler;

                size_t n = dummyDer_->reducedFun_->Domain(); // total variable count
                size_t m = dummyDer_->reducedFun_->Range(); // equation count

                std::vector<CG<Base> > indep(n);
                handler.makeVariables(indep);

                std::vector<CG<Base> > dep = dummyDer_->reducedFun_->Forward(0, indep);

                size_t algebraicCount = 0;
                for (size_t j = 0; j < dummyDer_->eqDifferentialInfo_.size(); j++) {
                    if (!dummyDer_->eqDifferentialInfo_[j]) {
                        algebraicCount++;
                    }
                }
                size_t stateCount = dummyDer_->eqDifferentialInfo_.size() - algebraicCount;

                /**
                 * Short independent variable vector (states will be considered constant)
                 */
                std::vector<AD<Base> > indepShort(n - stateCount);
                size_t pos = 0;
                for (size_t j = 0; j < dummyDer_->eqDifferentialInfo_.size(); j++) {
                    if (!dummyDer_->eqDifferentialInfo_[j]) {
                        indepShort[pos] = dummyDer_->x_[j];
                        pos++;
                    }
                }
                assert(pos == algebraicCount);
                for (size_t j = pos; j < dummyDer_->enodes_.size(); j++) {
                    indepShort[pos] = 0.0; // differential variable
                }
                Independent(indepShort);


                std::vector<AD<Base> > indep2(n);
                pos = 0;
                for (size_t j = 0; j < dummyDer_->eqDifferentialInfo_.size(); j++) {
                    if (!dummyDer_->eqDifferentialInfo_[j]) {
                        indep2[j] = indepShort[pos];
                        // algebraic variable normalization constant
                        normindep_[pos] = dummyDer_->norm_[j];
                        pos++;
                    } else {
                        indep2[j] = dummyDer_->x_[j]; // constant value
                    }
                }
                assert(pos == algebraicCount); // purely algebraic equations
                for (size_t j = pos; j < indepShort.size(); j++) {
                    indep2[j + stateCount] = indepShort[j]; // differential variable
                }

                // normalization constants for differential variables
                for (size_t j = 0; j < dummyDer_->vnodes_.size(); j++) {
                    if (dummyDer_->vnodes_[j]->derivativeOf() == NULL) {
                        Vnode<Base>* vDiff = dummyDer_->vnodes_[j]->derivative();
                        while (vDiff != NULL) {
                            normindep_[vDiff->index() - stateCount] = dummyDer_->norm_[j];
                            vDiff = vDiff->derivative();
                        }
                    }
                }

                Evaluator<Base, Base> evaluator(handler, dep);
                std::vector<AD<Base> > depNew = evaluator.evaluate(indep2);

                // turn every equation to a residual
                for (size_t i = 0; i < dummyDer_->eqDifferentialInfo_.size(); i++) {
                    if (dummyDer_->eqDifferentialInfo_[i]) {
                        Vnode<Base>* vDiff = dummyDer_->vnodes_[i]->derivative();
                        Enode<Base>* eq = dummyDer_->enodes_[i];

                        while (eq != NULL) {
                            assert(vDiff != NULL);
                            depNew[eq->index()] -= indepShort[vDiff->index() - stateCount];
                            normdep_[eq->index()] = dummyDer_->norm_[i];

                            eq = eq->derivative();
                            vDiff = vDiff->derivative();
                        }
                    }
                }

                reducedFunB_ = new ADFun<Base > (indepShort, depNew);
                assert(indepShort.size() == depNew.size());

                /**
                 * save new sparsity information
                 */
                jac_sparsity_ = jacobianReverseSparsitySet(*reducedFunB_);

                generateSparsityIndexes(jac_sparsity_, row_, col_);

                jac_.resize(row_.size());
            }

            ~Functor() {
                delete reducedFunB_;
            }

            int inputs() const {
                return reducedFunB_->Domain();
            }

            int values() const {
                return reducedFunB_->Range();
            }

            int operator()(const VectorB &x, VectorB & fvec) const {
                std::vector<Base> indep(x.rows()); //TODO: check this
                for (size_t j = 0; j < indep.size(); j++) {
                    indep[j] = x(j) * normindep_[j];
                }

                std::vector<Base> dep = reducedFunB_->Forward(0, indep);
                for (size_t j = 0; j < dep.size(); j++) {
                    fvec(j) = dep[j] / normdep_[j];
                }

                return 0;
            }

            int df(const VectorB &x, MatrixB & fjac) {

                std::vector<Base> indep(x.rows()); //TODO: check this
                for (size_t j = 0; j < indep.size(); j++) {
                    indep[j] = x(j) * normindep_[j];
                }

                size_t n_sweep = reducedFunB_->SparseJacobianReverse(indep, jac_sparsity_,
                                                                     row_, col_, jac_, work_);

                for (size_t pos = 0; pos < jac_.size(); pos++) {
                    size_t i = row_[pos];
                    size_t j = col_[pos];
                    fjac(i, j) = jac_[pos] / normdep_[j] * normindep_[i];
                }

                return 0;
            }
        };

        static Base evalCondition(const MatrixB& mat) {
            /**
             * determine eigenvalues
             */
            VectorCB eigenv = mat.eigenvalues();
#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "    eigen values:\n" << eigenv << "\n";
#endif

            /**
             * determine condition
             */
            if (eigenv(0).imag() != 0) {
                return std::numeric_limits<Base>::quiet_NaN();
            }
            Base max = std::abs(eigenv(0).real());
            Base min = max;

            for (size_t r = 1; r < eigenv.rows(); r++) {
                if (eigenv(r).imag() != 0) {
                    return std::numeric_limits<Base>::quiet_NaN();
                }
                Base eigv = std::abs(eigenv(r).real());
                if (eigv > max) {
                    max = eigv;
                } else if (eigv < min) {
                    min = eigv;
                }
            }

            // the condition number
            return max / min;
        }

    };
}


#endif

