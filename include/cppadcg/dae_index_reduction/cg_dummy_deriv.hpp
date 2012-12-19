#ifndef CPPAD_CG_DUMMY_DERIV_INCLUDED
#define	CPPAD_CG_DUMMY_DERIV_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <cppadcg/dae_index_reduction/cg_pantelides.hpp>

#include "cg_dae_var_info.hpp"

//#include <unsupported/Eigen/NonLinearOptimization>
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
        typedef CG<Base> CGBase;
        typedef AD<CGBase> ADCG;
        typedef Eigen::Matrix<Base, Eigen::Dynamic, 1 > VectorB;
        typedef Eigen::Matrix<std::complex<Base>, Eigen::Dynamic, 1 > VectorCB;
        typedef Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic> MatrixB;
    protected:
        // typical values;
        std::vector<Base> x_;
        // normalization constants for the variables (in the original order)
        std::vector<Base> normVar_;
        // normalization constants for the equations
        std::vector<Base> normEq_;
        /** Jacobian sparsity pattern of the reduced system
         * (in the original variable order)
         */
        std::vector<bool> jacSparsity_;
        // the initial index of time derivatives
        size_t diffVarStart_;
        // the initial index of the differentiated equations
        size_t diffEqStart_;
        /** normalized Jacobian of the index one system's  differentiated
         *  equations relative to the time derivatives
         *  (in the new variable order).
         */
        Eigen::SparseMatrix<Base, Eigen::RowMajor> jacobian_;
        /**
         * Dummy derivatives
         */
        std::set<Vnode<Base>* > dummyD_;
        /**
         * 
         */
        std::vector<DaeVarInfo> reducedVarInfo_;
    public:

        DummyDerivatives(ADFun<CG<Base> >* fun,
                         const std::vector<DaeVarInfo>& varInfo,
                         const std::vector<Base>& x,
                         const std::vector<Base>& normVar,
                         const std::vector<Base>& normEq) :
            Plantelides<Base>(fun, varInfo),
            x_(x),
            normVar_(normVar),
            normEq_(normEq),
            diffVarStart_(0),
            diffEqStart_(fun->Range()) {

            typename std::vector<Vnode<Base>*> ::const_iterator j;
            for (j = this->vnodes_.begin(); j != this->vnodes_.end(); ++j) {
                Vnode<Base>* jj = *j;
                if (jj->derivativeOf() != NULL) {
                    diffVarStart_ = jj->index();
                    break;
                }
            }
        }

        virtual inline std::vector<bool> timeDerivativeVariables() {
            std::vector<bool> tderiv(this->vnodes_.size());
            for (size_t j = 0; j< this->vnodes_.size(); j++) {
                Vnode<Base>* jj = this->vnodes_[j];
                bool isDeriv = jj->derivativeOf() != NULL && dummyD_.find(jj) == dummyD_.end();
                tderiv[jj->tapeIndex()] = isDeriv;
            }

            return tderiv;
        }

        virtual inline ADFun<CG<Base> >* reduceIndex(std::vector<DaeVarInfo>& newVarInfo) throw (CGException) {
            Plantelides<Base>::reduceIndex(reducedVarInfo_);

            if (this->reducedFun_ == NULL)
                return NULL; //nothing to do

            //solveDAESystem();

            determineJacobian();


            // variables of interest
            std::vector<Vnode<Base>* > vars;
            vars.reserve(this->vnodes_.size() - diffVarStart_);
            typename std::vector<Vnode<Base>* >::const_reverse_iterator rj;
            for (rj = this->vnodes_.rbegin(); rj != this->vnodes_.rend(); ++rj) {
                Vnode<Base>* jj = *rj;
                if (jj->derivativeOf() != NULL && jj->derivative() == NULL) {
                    vars.push_back(jj); // highest order time derivatives in the index 1 model
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


            Eigen::SparseMatrix<Base> workJac;

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

                // create the Jacobian for the selected variables and equations
                workJac.setZero();
                workJac.resize(eqs.size(), vars.size());
                for (size_t i = 0; i < eqs.size(); i++) {
                    Enode<Base>* ii = eqs[i];
                    for (size_t j = 0; j < vars.size(); j++) {
                        Vnode<Base>* jj = vars[j];
                        workJac.coeffRef(i, j) = jacobian_.coeff(ii->index() - diffEqStart_, jj->index() - diffVarStart_);
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
            }


            /**
             * Output information
             */
            typename std::set<Vnode<Base>* >::const_iterator j;
            for (j = dummyD_.begin(); j != dummyD_.end(); ++j) {
                assert((*j)->derivativeOf() != NULL);
                reducedVarInfo_[(*j)->derivativeOf()->tapeIndex()] = -1;
            }

            newVarInfo = reducedVarInfo_;

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "## dummy derivatives:\n";

            for (j = dummyD_.begin(); j != dummyD_.end(); ++j)
                std::cout << "# " << **j << "\n";
            std::cout << "# \n";
#endif

            return this->reducedFun_;
        }

        inline ADFun<CG<Base> >* reduceEquations(std::vector<DaeVarInfo>& newVarInfo) {
            using namespace std;
            using std::vector;

            if (this->reducedFun_ == NULL) {
                return NULL; // nothing to do
            }

            /**
             * Generate an operation graph
             */
            CodeHandler<Base> handler;

            vector<CGBase> indep0(this->reducedFun_->Domain());
            handler.makeVariables(indep0);

            const vector<CGBase> res0 = this->reducedFun_->Forward(0, indep0);

            /**
             * determine if the time derivative variable is present in the 
             * equation only multiplied by a constant value
             */
            set<size_t> removedEquations;
            set<Vnode<Base>* > removedDummy;
            typename set<Vnode<Base>* >::const_iterator j;
            for (j = dummyD_.begin(); j != dummyD_.end(); ++j) {
                Vnode<Base>* dummy = *j;
                Enode<Base>* i = dummy->assigmentEquation();

                try {
                    // eliminate all references to the dummy variable by substitution
                    handler.substituteIndependent(indep0[dummy->tapeIndex()], res0[i->index()]);

                    removedEquations.insert(i->index());
                    removedDummy.insert(dummy);
                } catch (const CGException& ex) {
                    // unable to solve for a dummy variable: keep the equation and variable
                }
            }

            vector<size_t> orig2New(this->vnodes_.size(), -1);
            vector<bool> usedOrig(this->vnodes_.size(), false);
            vector<Vnode<Base>* > vnodesNew;
            typename vector<Vnode<Base>* >::const_iterator j2;
            for (j2 = this->vnodes_.begin(); j2 != this->vnodes_.end(); ++j2) {
                if (removedDummy.find(*j2) == removedDummy.end()) {
                    vnodesNew.push_back(*j2);
                    orig2New[(*j2)->tapeIndex()] = (*j2)->index();
                    usedOrig[(*j2)->tapeIndex()] = true;
                }
            }

            // determine the new indexes in the tape
            for (size_t p = 0; p < orig2New.size(); p++) {
                if (!usedOrig[p]) {
                    for (size_t p2 = p + 1; p2 < orig2New.size(); p2++) {
                        if (usedOrig[p2]) {
                            orig2New[p2]--;
                        }
                    }
                }
            }

            /**
             * create a new tape without the dummy derivatives and with 
             * less equations
             */
            // less variables
            vector<ADCG> indepNew(handler.getIndependentVariableSize());
            assert(indepNew.size() == vnodesNew.size());

            for (size_t p = 0; p < orig2New.size(); p++) {
                if (usedOrig[p]) {
                    indepNew[orig2New[p]] = x_[p];
                }
            }

            Independent(indepNew);

            // less equations
            vector<CGBase> resNew(this->reducedFun_->Range() - removedDummy.size());
            for (size_t i = 0, p = 0; i< this->enodes_.size(); i++) {
                if (removedEquations.find(i) == removedEquations.end()) {
                    resNew[p] = res0[i];
                    p++;
                }
            }

            Evaluator<Base, CGBase> evaluator0(handler, resNew);
            vector<ADCG> depNew = evaluator0.evaluate(indepNew);

            ADFun<CG<Base> >* shortFun = new ADFun<CGBase > (indepNew, depNew);
            try {
                vector<string> vnames(vnodesNew.size());
                for (size_t p = 0; p < vnodesNew.size(); p++) {
                    Vnode<Base>* n = vnodesNew[p];
                    vnames[orig2New[n->tapeIndex()]] = n->name();
                }

#ifdef CPPAD_CG_DAE_VERBOSE
                this->printModel(shortFun, vnames);
#endif

                /**
                 * Output information
                 */
                assert(usedOrig.size() == reducedVarInfo_.size());

                newVarInfo = reducedVarInfo_;
                for (int p = usedOrig.size() - 1; p >= 0; p--) {
                    if (!usedOrig[p]) {
                        newVarInfo.erase(newVarInfo.begin() + p);
                        for (size_t pp = 0; pp < usedOrig.size(); pp++) {
                            if (newVarInfo[pp].getDerivativeOf() > p) {
                                newVarInfo[pp] = DaeVarInfo(newVarInfo[pp].getDerivativeOf() - 1);
                            }
                        }
                    }
                }

            } catch (...) {
                delete shortFun;
                throw;
            }

            return shortFun;
        }

        inline virtual ~DummyDerivatives() {
        }

    protected:

        inline void solveDAESystem() {
            throw 1; // not finished!!!!
            /**
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
             **/
        }

        /**
         * Determines the Jacobian relative to the differential variables
         * (e.g. dxdt)
         */
        inline void determineJacobian() {
            using namespace std;
            using std::vector;

            const size_t n = this->reducedFun_->Domain();
            const size_t m = this->reducedFun_->Range();

            jacSparsity_ = jacobianReverseSparsity < vector<bool>, CGBase > (*this->reducedFun_); // in the original variable order

            vector<size_t> row, col;
            row.reserve((this->vnodes_.size() - diffVarStart_) * (m - diffEqStart_));
            col.reserve(row.capacity());

            for (size_t i = diffEqStart_; i < m; i++) {
                for (size_t j = diffVarStart_; j < n; j++) {
                    assert(this->vnodes_[j]->derivativeOf() != NULL);
                    size_t t = this->vnodes_[j]->tapeIndex();
                    if (jacSparsity_[i * n + t]) {
                        row.push_back(i);
                        col.push_back(t);
                    }
                }
            }

            vector<CG<Base> > jac(row.size());

            vector<CG<Base> > indep(n);
            std::copy(x_.begin(), x_.end(), indep.begin());
            std::fill(indep.begin() + x_.size(), indep.end(), 0);

            CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
            this->reducedFun_->SparseJacobianReverse(indep, jacSparsity_,
                                                     row, col, jac, work);

            // resize and zero matrix
            jacobian_.resize(m - diffEqStart_, n - diffVarStart_);

            map<size_t, Vnode<Base>*> origIndex2var;
            for (size_t j = diffVarStart_; j< this->vnodes_.size(); j++) {
                Vnode<Base>* jj = this->vnodes_[j];
                origIndex2var[jj->tapeIndex()] = jj;
            }

            // normalize values
            for (size_t e = 0; e < jac.size(); e++) {
                Enode<Base>* eqOrig = this->enodes_[row[e]]->originalEquation();
                Vnode<Base>* vOrig = origIndex2var[col[e]]->originalVariable(this->fun_->Domain());

                // normalized jacobian value
                Base normVal = jac[e].getParameterValue() * normVar_[vOrig->tapeIndex()]
                        / normEq_[eqOrig->index()];

                size_t i = row[e]; // same order
                size_t j = origIndex2var[col[e]]->index(); // different order than in model/tape

                jacobian_.coeffRef(i - diffEqStart_, j - diffVarStart_) = normVal;
            }

            jacobian_.makeCompressed();

#ifdef CPPAD_CG_DAE_VERBOSE
            cout << "partial jacobian:\n" << jacobian_ << "\n\n";
            //cout << jacobian_.triangularView<Eigen::Lower > () << "\n\n";
#endif
        }

        inline void selectDummyDerivatives(const std::vector<Enode<Base>* >& eqs,
                                           const std::vector<Vnode<Base>* >& vars,
                                           Eigen::SparseMatrix<Base>& subsetJac) throw (CGException) {

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

            /**
             * Fill in the Jacobian subset for the selected equations and variables
             */
            subsetJac.resize(eqs.size(), vars.size());
            std::vector<size_t> rowNnz(eqs.size()); // the number of non-zero elements per row
            std::vector<size_t> rowNnzCol(eqs.size()); // the last defined column for each row
            for (size_t i = 0; i < eqs.size(); i++) {
                Enode<Base>* ii = eqs[i];
                for (size_t j = 0; j < vars.size(); j++) {
                    Vnode<Base>* jj = vars[j];
                    Base val = jacobian_.coeff(ii->index() - diffEqStart_, jj->index() - diffVarStart_);
                    if (val != Base(0.0)) {
                        subsetJac.coeffRef(i, j) = val;
                        rowNnz[i]++;
                        rowNnzCol[i] = j;
                    }
                }
            }
#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "subset Jac:\n" << subsetJac << "\n";
#endif

            MatrixB workJac(eqs.size(), eqs.size());

            /**
             * Determine the columns that cannot be removed
             */
            std::set<size_t> fixedCols;
            for (size_t i = 0; i < rowNnz.size(); ++i) {
                if (rowNnz[i] == 1) {
                    fixedCols.insert(rowNnzCol[i]);
                }
            }

#ifdef CPPAD_CG_DAE_VERBOSE
            if (!fixedCols.empty()) {
                std::cout << " fixed columns:";
                for (std::set<size_t>::const_iterator it = fixedCols.begin(); it != fixedCols.end(); ++it) {
                    std::cout << " " << *vars[*it];
                }
                std::cout << "\n";
            }
#endif

            /**
             * column indexes that can be added/removed from the selection
             */
            std::vector<size_t> freeCols;
            for (size_t j = 0; j < vars.size(); ++j) {
                if (fixedCols.find(j) == fixedCols.end()) {
                    freeCols.push_back(j);
                }
            }

            std::vector<size_t> vcols2keep(eqs.size() - fixedCols.size());
            for (size_t c = 0; c < vcols2keep.size(); c++) {
                vcols2keep[c] = c;
            }

            // number of columns/variables to remove (the remaining will be dummy derivatives)
            std::vector<size_t> cols2keep(eqs.size());
            {
                std::set<size_t> cols2keepAux(fixedCols);
                for (size_t c = 0; c < vcols2keep.size(); c++) {
                    cols2keepAux.insert(freeCols[c]);
                }
                std::copy(cols2keepAux.begin(), cols2keepAux.end(), cols2keep.begin());
            }

            /**
             * Brute force approach!!!
             */
            std::vector<size_t> bestCols2keep;
            Base bestCond = std::numeric_limits<Base>::max();
            size_t bestTotalOrder = 0;

            while (true) {

#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << " ## column selection: ";
                for (size_t s = 0; s < cols2keep.size(); s++)
                    std::cout << cols2keep[s] << " ";
                std::cout << " \n";
#endif
                workJac.setZero(eqs.size(), eqs.size());
                for (size_t c = 0; c < cols2keep.size(); ++c) {
                    typename Eigen::SparseMatrix<Base>::InnerIterator itCol(subsetJac, cols2keep[c]);
                    for (; itCol; ++itCol) {
                        assert(itCol.col() == int(cols2keep[c]));
                        workJac(itCol.row(), c) = itCol.value();
                    }
                }

#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << "    current jac:\n" << workJac << "\n";
#endif

                Base cond = evalBestMatrixCondition(workJac);

#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << "    condition: " << cond << "\n";
#endif

                if (cond == cond) {
                    // not NaN
                    size_t totalOrd = 0;
                    for (size_t j = 0; j < cols2keep.size(); j++) {
                        totalOrd += vars[cols2keep[j]]->order();
                    }
                    if ((totalOrd > bestTotalOrder && cond / Base(10.0) <= bestCond) ||
                            (totalOrd == bestTotalOrder && cond < bestCond) ||
                            (totalOrd < bestTotalOrder && cond * Base(10.0) <= bestCond)) {
                        bestTotalOrder = totalOrd;
                        bestCond = cond;
                        bestCols2keep = cols2keep;
                    }
                }

                /**
                 * determine the next set of columns
                 */
                cols2keep = nextColumnSelection(fixedCols, freeCols, vcols2keep);
                if (cols2keep.empty())
                    break;
            };

            if (bestCols2keep.empty()) {
                throw CGException("Failed to select dummy derivatives! "
                                  "The resulting system is probably singular for the provided data.");
            }

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "## new dummy derivatives (condition = " << bestCond << "): ";
            for (size_t c = 0; c < bestCols2keep.size(); c++)
                std::cout << *vars[bestCols2keep[c]] << "; ";
            std::cout << " \n\n";
#endif

            for (size_t c = 0; c < bestCols2keep.size(); c++) {
                dummyD_.insert(vars[bestCols2keep[c]]);
            }

        }

        /**
         * 
         * \param fixedCols Column indeces that must be selected
         * \param freeCols Column that can be selected (excluding the fixedCols)
         * \param vcols2keep The previous column selection from the free columns
         * @return the next column selection
         */
        inline std::vector<size_t > nextColumnSelection(const std::set<size_t>& fixedCols,
                                                        const std::vector<size_t>& freeCols,
                                                        std::vector<size_t>& vcols2keep) const {

            if (vcols2keep.empty()) {
                return std::vector<size_t > (0); // end of combinations
            }

            if (vcols2keep.back() == freeCols.size() - 1) {
                if (vcols2keep[0] == freeCols.size() - vcols2keep.size())
                    return std::vector<size_t > (0); // end of combinations

                for (size_t cc = 1; cc < vcols2keep.size(); cc++) {
                    if (vcols2keep[cc] == freeCols.size() - (vcols2keep.size() - cc)) {
                        vcols2keep[cc - 1]++;
                        for (size_t cc2 = cc; cc2 < vcols2keep.size(); cc2++) {
                            vcols2keep[cc2] = vcols2keep[cc2 - 1] + 1;
                        }
                        break;
                    }
                }
            } else {
                vcols2keep.back()++;
            }

            std::set<size_t> cols2keep(fixedCols);

            for (size_t c = 0; c < vcols2keep.size(); c++) {
                size_t vColIndex = freeCols[vcols2keep[c]];
                cols2keep.insert(vColIndex);
            }

            return std::vector<size_t > (cols2keep.begin(), cols2keep.end());
        }

        inline void generateSemiExplicitDAE() {
            /**
             * 
             * TODO!!!
             * 
             */
            using std::vector;

            CodeHandler<Base> handler;

            vector<CGBase> indep0(this->reducedFun_->Domain());
            handler.makeVariables(indep0);

            const vector<CGBase> dep0 = this->reducedFun_->Forward(0, indep0);

            const vector<bool> algEq(dep0.size());

            /**
             * make relations between dependent and independent variables.
             * Equations may be removed when used to calculate dummy derivatives.
             */

            // check if it is possible to get a semi-explicit DAE


            for (size_t diffj = diffVarStart_; diffj < this->vnodes_.size(); diffj++) {
                // find equation used to determine it
                Vnode<Base>* diffjj = this->vnodes_[diffj];

                //   get original var -> get diff equation -> get norder equation
                Vnode<Base>* origj = diffjj->originalVariable();

                Enode<Base>* eq = this->enodes_[origj->index()];
                if (eq->isAlgebraic()) {
                    std::stringstream ss;
                    ss << "Unable to produce a semi-explicit DAE system due to the presence"
                            " of the algebraic variable '" << origj->index() << "' in new equation(s)"
                            " generated by differentiation of existing algebraic equations.";
                    throw CGException(ss.str());
                }

                size_t order = diffjj->order();
                for (size_t o = 0; o < order - 1; o++) {
                    eq = eq->derivative();
                }

                SourceCodeFragment<Base>* code = dep0[eq->index()].getSourceCodeFragment();

            }

#ifdef CPPAD_CG_DAE_VERBOSE
            this->printModel(this->reducedFun_);
#endif
        }

        /**
         * Determines the best matrix
         * \param mat The matrix
         * @return The best condition value (lowest possible and real)
         */
        static Base evalBestMatrixCondition(const MatrixB& mat) {

            Eigen::FullPivLU<MatrixB> lu = mat.fullPivLu();
            //  MatrixB l = MatrixB::Identity(mat.rows(), mat.cols());
            //  l.template triangularView<Eigen::StrictlyLower > () = lu.matrixLU();
            MatrixB u = lu.matrixLU().template triangularView<Eigen::Upper > ();

            //  std::cout << "mat:\n" << mat << "\n\n";
            //  std::cout << "L:\n" << l << "\n\n";
            //  std::cout << "U:\n" << u << "\n\n";

            //VectorCB eigenv = u.eigenvalues();            
            //std::cout << "    eigen values:\n" << eigenv << "\n";

            /**
             * determine condition of U 
             * (the eigenvalues are in the diagonal)
             */
            if (u(0, 0) == 0) {
                return std::numeric_limits<Base>::quiet_NaN();
            }
            Base max = std::abs(u(0, 0));
            Base min = max;

            for (int r = 1; r < u.rows(); r++) {
                if (u(r, r) == 0) {
                    return std::numeric_limits<Base>::quiet_NaN();
                }
                Base eigv = std::abs(u(r, r));
                if (eigv > max) {
                    max = eigv;
                } else if (eigv < min) {
                    min = eigv;
                }
            }

            // the condition number
            return max / min;
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
                        normindep_[pos] = dummyDer_->normVar_[j];
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
                            normindep_[vDiff->index() - stateCount] = dummyDer_->normVar_[j];
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
                            normdep_[eq->index()] = dummyDer_->normVar_[i];

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

    };
}

#endif
