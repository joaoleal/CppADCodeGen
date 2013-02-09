#ifndef CPPAD_CG_DUMMY_DERIV_INCLUDED
#define CPPAD_CG_DUMMY_DERIV_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <cppadcg/dae_index_reduction/cg_pantelides.hpp>
//#include <unsupported/Eigen/NonLinearOptimization>

namespace CppAD {

    /**
     * Sorts variable nodes according to the variable differentiation order
     * 
     * \param i
     * \param j
     * \return true if i should come before j
     */
    template<class Base>
    bool sortVnodesByOrder(Vnode<Base>* i, Vnode<Base>* j) {
        return (i->order() > j->order());
    }

    /**
     * Utility class used to sort varibles in the DAE
     */
    class DAEVarOrderInfo {
    public:
        size_t originalIndex;
        size_t originalIndex0;
        bool hasDerivatives;
        int order;

        inline DAEVarOrderInfo() :
            originalIndex(0),
            originalIndex0(0),
            hasDerivatives(false),
            order(-1) {
        }

        inline DAEVarOrderInfo(size_t moriginalIndex, size_t moriginalIndex0, bool mhasDerivatives, int morder) :
            originalIndex(moriginalIndex),
            originalIndex0(moriginalIndex0),
            hasDerivatives(mhasDerivatives),
            order(morder) {
        }
    };

    /**
     * Sorts variables based on the differentiation order, whether they are 
     * algebraic or differential and the order in the original model
     * 
     * \param i
     * \param j
     * \return true if i should come before j
     */
    inline bool sortVariablesByOrder(const DAEVarOrderInfo& i, const DAEVarOrderInfo& j) {
        if (j.order < i.order) {
            return true;
        } else if (j.order > i.order) {
            return false;
        } else if (i.hasDerivatives == j.hasDerivatives) {
            return j.originalIndex > i.originalIndex;
        } else {
            return i.hasDerivatives;
        }
    }

    /**
     * Sorts equations according to the index/order of the assigned variables
     * 
     * \param i
     * \param j
     * \return true if i should come before j
     */
    inline bool sortEquationByAssignedOrder(const DaeEquationInfo& i, const DaeEquationInfo& j) {
        if (i.getAssignedVarIndex() < 0) {
            if (j.getAssignedVarIndex() < 0)
                return i.getOriginalIndex() < j.getOriginalIndex();
            else
                return false;
        } else {
            if (j.getAssignedVarIndex() >= 0)
                return i.getAssignedVarIndex() < j.getAssignedVarIndex();
            else
                return true;
        }
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
        // normalization constants for the variables (in the original order)
        std::vector<Base> normVar_;
        // normalization constants for the equations
        std::vector<Base> normEq_;
        /** 
         * Jacobian sparsity pattern of the reduced system
         * (in the original variable order)
         */
        std::vector<bool> jacSparsity_;
        // the initial index of time derivatives
        size_t diffVarStart_;
        // the initial index of the differentiated equations
        size_t diffEqStart_;
        /**
         * Normalized Jacobian of the index one system's  differentiated
         * equations relative to the time derivatives
         * (in the new variable order).
         */
        Eigen::SparseMatrix<Base, Eigen::RowMajor> jacobian_;
        /**
         * Dummy derivatives
         */
        std::set<Vnode<Base>* > dummyD_;
        bool reduceEquations_;
        bool generateSemiExplicitDae_;
        bool reorder_;
    public:

        /**
         * Creates the DAE index reduction algorithm that implements the dummy
         * derivatives method.
         * 
         * \param fun The DAE model
         * \param varInfo DAE model variable classification
         * \param x typical variable values (used to determine Jacobian values)
         * \param normVar variable normalization values
         * \param normEq equation normalization values
         */
        DummyDerivatives(ADFun<CG<Base> >* fun,
                         const std::vector<DaeVarInfo>& varInfo,
                         const std::vector<Base>& x,
                         const std::vector<Base>& normVar,
                         const std::vector<Base>& normEq) :
            Plantelides<Base>(fun, varInfo, x),
            normVar_(normVar),
            normEq_(normEq),
            diffVarStart_(0),
            diffEqStart_(fun->Range()),
            reduceEquations_(true),
            generateSemiExplicitDae_(false),
            reorder_(true) {

            typename std::vector<Vnode<Base>*> ::const_iterator j;
            for (j = this->vnodes_.begin(); j != this->vnodes_.end(); ++j) {
                Vnode<Base>* jj = *j;
                if (jj->antiDerivative() != NULL) {
                    diffVarStart_ = jj->index();
                    break;
                }
            }
        }

        inline bool isGenerateSemiExplicitDae() const {
            return generateSemiExplicitDae_;
        }

        inline void setGenerateSemiExplicitDae(bool generateSemiExplicitDae) {
            generateSemiExplicitDae_ = generateSemiExplicitDae;
        }

        inline bool isReduceEquations() const {
            return reduceEquations_;
        }

        inline void setReduceEquations(bool reduceEquations) {
            reduceEquations_ = reduceEquations;
        }

        inline bool isReorder() const {
            return reorder_;
        }

        inline void setReorder(bool reorder) {
            reorder_ = reorder;
        }

        /**
         * Performs the DAE differentiation index reductions
         * 
         * \param newVarInfo Variable related information of the reduced index
         *                   model
         * \param equationInfo Equation related information of the reduced index
         *                     model
         * \return the reduced index model (must be deleted by user)
         */
        virtual inline ADFun<CG<Base> >* reduceIndex(std::vector<DaeVarInfo>& newVarInfo,
                                                     std::vector<DaeEquationInfo>& newEqInfo) throw (CGException) {

            /**
             * Variable information for the reduced
             */
            std::vector<DaeVarInfo> reducedVarInfo;
            /**
             * Equation information for the reduced model
             */
            std::vector<DaeEquationInfo> reducedEqInfo;

            std::auto_ptr<ADFun<CG<Base> > > fun(Plantelides<Base>::reduceIndex(reducedVarInfo, reducedEqInfo));
            if (fun.get() == NULL)
                return NULL; //nothing to do (no index reduction required)

            newEqInfo = reducedEqInfo; // copy
            addDummyDerivatives(reducedVarInfo, newVarInfo);

            if (reduceEquations_) {
                std::auto_ptr<ADFun<CG<Base> > > funShort = reduceEquations(reducedVarInfo, newVarInfo,
                                                                            reducedEqInfo, newEqInfo);
                fun = funShort;
            }

            if (generateSemiExplicitDae_) {
                std::vector<DaeVarInfo> varInfo = newVarInfo; // copy
                std::vector<DaeEquationInfo> eqInfo = newEqInfo; // copy
                std::auto_ptr<ADFun<CG<Base> > > semiExplicit = generateSemiExplicitDAE(fun.get(),
                                                                                        varInfo, newVarInfo,
                                                                                        eqInfo, newEqInfo);
                fun = semiExplicit;
            }

            if (reorder_) {
                std::vector<DaeVarInfo> varInfo = newVarInfo; // copy
                std::vector<DaeEquationInfo> eqInfo = newEqInfo; // copy
                std::auto_ptr<ADFun<CG<Base> > > reorderedFun = reorderModelEqNVars(fun.get(),
                                                                                    varInfo, newVarInfo,
                                                                                    eqInfo, newEqInfo);
                fun = reorderedFun;
            }

            return fun.release();
        }

        inline virtual ~DummyDerivatives() {
        }

    protected:

        virtual inline void addDummyDerivatives(const std::vector<DaeVarInfo>& varInfo,
                                                std::vector<DaeVarInfo>& newVarInfo) throw (CGException) {

            determineJacobian();

            // variables of interest
            std::vector<Vnode<Base>* > vars;
            vars.reserve(this->vnodes_.size() - diffVarStart_);
            typename std::vector<Vnode<Base>* >::const_reverse_iterator rj;
            for (rj = this->vnodes_.rbegin(); rj != this->vnodes_.rend(); ++rj) {
                Vnode<Base>* jj = *rj;
                if (jj->antiDerivative() != NULL && jj->derivative() == NULL) {
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
                    Vnode<Base>* v = (*j)->antiDerivative();
                    if (v != NULL && v->antiDerivative() != NULL) {
                        varsNew.push_back(v);
                    }
                }
                vars.swap(varsNew);
            }


            /**
             * Prepare the output information
             */
            newVarInfo = varInfo; //copy
            typename std::set<Vnode<Base>* >::const_iterator j;
            for (j = dummyD_.begin(); j != dummyD_.end(); ++j) {
                assert((*j)->antiDerivative() != NULL);
                assert((*j)->tapeIndex() >= 0);
                newVarInfo[(*j)->tapeIndex()].setAntiDerivative(-1);
            }

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "## dummy derivatives:\n";

            for (j = dummyD_.begin(); j != dummyD_.end(); ++j)
                std::cout << "# " << **j << "   " << newVarInfo[(*j)->tapeIndex()].getName() << "\n";
            std::cout << "# \n";
            Plantelides<Base>::printModel(this->reducedFun_, newVarInfo);
#endif

        }

        /**
         * Attempts to reduce the number of equations by variable substitution.
         * 
         * \param newVarInfo Variable information of the resulting model
         * \return The new DAE reduced model with (possibly) less equations and
         *         variables
         */
        inline std::auto_ptr<ADFun<CGBase > > reduceEquations(const std::vector<DaeVarInfo>& reducedVarInfo,
                                                              std::vector<DaeVarInfo>& newVarInfo,
                                                              const std::vector<DaeEquationInfo>& reducedEqInfo,
                                                              std::vector<DaeEquationInfo>& newEqInfo) {
            using namespace std;
            using std::vector;

            /**
             * Generate an operation graph
             */
            CodeHandler<Base> handler;

            vector<CGBase> indep0(this->reducedFun_->Domain());
            handler.makeVariables(indep0);

            vector<CGBase> res0 = this->reducedFun_->Forward(0, indep0);

            /**
             * maps the equations indexes of the reduced model to the new 
             * equation indexes in the model with less equations and variables
             * (removed equations have negative indexes)
             */
            std::vector<int> eqIndexReduced2Short(this->enodes_.size());
            for (size_t i = 0; i < eqIndexReduced2Short.size(); i++) {
                eqIndexReduced2Short[i] = i;
            }

            /**
             * maps the variables indexes in the tape of the reduced model to 
             * the  new tape indexes in the model with less equations and
             * variables (removed variables have negative indexes)
             */
            std::vector<int> tapeIndexReduced2Short(reducedVarInfo.size());
            for (size_t j = 0; j < tapeIndexReduced2Short.size(); j++) {
                tapeIndexReduced2Short[j] = j;
            }

            /**
             * attempt to eliminate dummy derivatives and the equation they are
             * assigned to
             */
#ifdef CPPAD_CG_DAE_VERBOSE
            std::set<size_t> erasedVariables;
            std::set<size_t> erasedEquations;
#endif

            typename set<Vnode<Base>* >::const_iterator j;
            for (j = dummyD_.begin(); j != dummyD_.end(); ++j) {
                Vnode<Base>* dummy = *j;
                Enode<Base>* i = dummy->assigmentEquation();

                try {
                    // eliminate all references to the dummy variable by substitution
                    handler.substituteIndependent(indep0[dummy->tapeIndex()], res0[i->index()]);
                    tapeIndexReduced2Short[dummy->tapeIndex()] = -1;
                    eqIndexReduced2Short[i->index()] = -1;

#ifdef CPPAD_CG_DAE_VERBOSE
                    std::cout << "######### use equation " << i->index() << " to solve for variable " << dummy->name() << std::endl;
                    erasedVariables.insert(dummy->tapeIndex());
                    erasedEquations.insert(i->index());
                    printModel(handler, res0, reducedVarInfo, erasedVariables, erasedEquations);
#endif
                } catch (const CGException& ex) {
                    // unable to solve for a dummy variable: keep the equation and variable
                }
            }

            // determine the new equation indexes
            for (size_t i = 0; i < eqIndexReduced2Short.size(); i++) {
                if (eqIndexReduced2Short[i] < 0) { // removed equation
                    for (size_t ii = i + 1; ii < eqIndexReduced2Short.size(); ii++) {
                        if (eqIndexReduced2Short[ii] >= 0) {
                            eqIndexReduced2Short[ii]--;
                        }
                    }
                }
            }

            // determine the new indexes in the tape
            for (size_t p = 0; p < tapeIndexReduced2Short.size(); p++) {
                if (tapeIndexReduced2Short[p] < 0) {
                    // removed from model
                    for (size_t p2 = p + 1; p2 < tapeIndexReduced2Short.size(); p2++) {
                        if (tapeIndexReduced2Short[p2] >= 0) {
                            tapeIndexReduced2Short[p2]--;
                        }
                    }
                }
            }

            /**
             * Prepare the output information
             */
            assert(tapeIndexReduced2Short.size() == reducedVarInfo.size());

            newVarInfo = reducedVarInfo; // copy
            for (int p = tapeIndexReduced2Short.size() - 1; p >= 0; p--) {
                if (tapeIndexReduced2Short[p] < 0) { // removed from model
                    newVarInfo.erase(newVarInfo.begin() + p);
                    for (size_t pp = 0; pp < tapeIndexReduced2Short.size(); pp++) {
                        DaeVarInfo& v = newVarInfo[pp];
                        if (v.getAntiDerivative() > p) {
                            v.setAntiDerivative(v.getAntiDerivative() - 1);
                        } else if (v.getAntiDerivative() == p) {
                            v.setAntiDerivative(-1);
                        }
                        if (v.getDerivative() > p) {
                            v.setDerivative(v.getDerivative() - 1);
                        } else if (v.getDerivative() == p) {
                            v.setDerivative(-1);
                        }
                    }
                }
            }

            newEqInfo = reducedEqInfo; // copy

            for (int p = eqIndexReduced2Short.size() - 1; p >= 0; p--) {
                if (eqIndexReduced2Short[p] < 0) {// removed from model
                    newEqInfo.erase(newEqInfo.begin() + p);
                } else {
                    DaeEquationInfo& eq = newEqInfo[p];
                    int reducedVIndex = eq.getAssignedVarIndex();
                    if (reducedVIndex >= 0)
                        eq.setAssignedVarIndex(tapeIndexReduced2Short[reducedVIndex]);
                    if (eq.getAntiDerivative() >= 0)
                        eq.setAntiDerivative(eqIndexReduced2Short[eq.getAntiDerivative()]);
                }
            }

            /**
             * Implement the model after after the reduction of equations and 
             * variables by substitution
             */
            std::auto_ptr<ADFun<CGBase > > shortFun(generateReorderedModel(handler, res0,
                                                                           reducedVarInfo, newVarInfo,
                                                                           reducedEqInfo, newEqInfo));

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "DAE with less equations and variables:\n";
            Plantelides<Base>::printModel(shortFun.get(), newVarInfo);
#endif

            return shortFun;
        }

        /**
         * Attempts to generate a semi-explicit DAE.
         * 
         * \param reorder place all the differential equations and variables
         *                together
         * \param differentialEqs 
         * \return The new semi-explicit DAE model with less variables (without
         *         the time derivative variables)
         */
        inline std::auto_ptr<ADFun<CGBase > > generateSemiExplicitDAE(ADFun<CG<Base> >* fun,
                                                                      const std::vector<DaeVarInfo>& varInfo,
                                                                      std::vector<DaeVarInfo>& newVarInfo,
                                                                      const std::vector<DaeEquationInfo>& eqInfo,
                                                                      std::vector<DaeEquationInfo>& newEqInfo) throw (CGException) {
            using namespace std;
            using std::vector;

            assert(fun != NULL);

            newEqInfo = eqInfo; // copy (we will have the same number of equations)

            /**
             * Generate an operation graph
             */
            CodeHandler<Base> handler;

            vector<CGBase> indep0(fun->Domain());
            handler.makeVariables(indep0);

            vector<CGBase> res0 = fun->Forward(0, indep0);

            /**
             * determine the variable indexes after the elimination of the time
             * derivatives
             */
            vector<int> varIndexOld2New(varInfo.size(), -1);
            size_t count = 0;
            for (size_t j = 0; j != varInfo.size(); ++j) {
                // exclude derivatives (they will be removed)
                if (varInfo[j].getAntiDerivative() < 0) {
                    varIndexOld2New[j] = count++;
                }
            }

            /**
             * Eliminate time derivatives from equations
             */
            for (size_t i = 0; i < eqInfo.size(); ++i) {
                const DaeEquationInfo& ii = eqInfo[i];
                int j = ii.getAssignedVarIndex();
                if (j < 0)
                    continue;
                const DaeVarInfo& jj = varInfo[j];

                if (jj.getAntiDerivative() >= 0) {
                    try {
                        CGBase& dep = res0[i]; // the equation residual
                        CGBase& indep = indep0[j]; // the time derivative

                        handler.substituteIndependent(indep, dep); // removes indep from the list of variables

                        SourceCodeFragment<Base>* alias = indep.getSourceCodeFragment();
                        assert(alias != NULL && alias->operation() == CGAliasOp);
                        dep.getSourceCodeFragment()->makeAlias(alias->arguments()[0]);

                        // it is now an explicit differential equation
                        newEqInfo[i].setExplicit(true);
                        // the derivative variable will disappear, associate the equation with the original variable
                        newEqInfo[i].setAssignedVarIndex(varIndexOld2New[jj.getAntiDerivative()]);
                    } catch (const CGException& ex) {
                        // unable to solve for a dummy variable: keep the equation and variable
                        throw CGException(string("Failed to generate semi-explicit DAE: ") + ex.what());
                    }
                }
            }

            /**
             * Prepare the output information
             */
            newVarInfo = varInfo;
            for (int j = newVarInfo.size() - 1; j >= 0; --j) {
                if (newVarInfo[j].getAntiDerivative() >= 0) {
                    // a derivative
                    newVarInfo.erase(newVarInfo.begin() + j);
                }
            }
            for (size_t j = 0; j < newVarInfo.size(); j++) {
                newVarInfo[j].setDerivative(-1); // no derivatives in tape
            }

            /**
             * Implement the reordering and derivative variable elimination in
             * the model
             */
            std::auto_ptr<ADFun<CGBase > > semiExplicitFun(generateReorderedModel(handler, res0, varInfo, newVarInfo, eqInfo, newEqInfo));

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "Semi-Eplicit DAE:\n";
            Plantelides<Base>::printModel(semiExplicitFun.get(), newVarInfo);
#endif

            return semiExplicitFun;
        }

        inline std::auto_ptr<ADFun<CGBase > > reorderModelEqNVars(ADFun<CG<Base> >* fun,
                                                                  const std::vector<DaeVarInfo>& varInfo,
                                                                  std::vector<DaeVarInfo>& newVarInfo,
                                                                  const std::vector<DaeEquationInfo>& eqInfo,
                                                                  std::vector<DaeEquationInfo>& newEqInfo) {

            using namespace std;
            using std::vector;

            assert(fun != NULL);

            /**
             * Determine the variables that have derivatives in the model
             */
            std::set<size_t> oldVarWithDerivatives; // indexes of old variables (before reordering) with derivatives
            for (size_t i = 0; i < eqInfo.size(); i++) {
                if (eqInfo[i].isExplicit() && eqInfo[i].getAssignedVarIndex() >= 0) {
                    oldVarWithDerivatives.insert(eqInfo[i].getAssignedVarIndex());
                }
            }

            if (oldVarWithDerivatives.empty()) {
                // no semi-explicit model generated
                for (size_t j = 0; j < varInfo.size(); j++) {
                    int index = j;
                    bool differential = false;
                    while (varInfo[index].getAntiDerivative() >= 0) {
                        index = varInfo[index].getAntiDerivative();
                        differential = true;
                    }

                    if (differential) {
                        oldVarWithDerivatives.insert(index);
                    }
                }
            }

            /**
             * sort variables
             */
            std::vector<DAEVarOrderInfo> varOrder(varInfo.size());
            for (size_t j = 0; j < varInfo.size(); j++) {
                size_t j0;
                int derivOrder = this->determineVariableDiffOrder(varInfo, j, j0);
                if (varInfo[j].isIntegratedVariable()) {
                    derivOrder = -2; // so that it goes last
                }
                bool hasDerivatives = oldVarWithDerivatives.find(j) != oldVarWithDerivatives.end();
                varOrder[j] = DAEVarOrderInfo(j, j0, hasDerivatives, derivOrder);
            }

            std::sort(varOrder.begin(), varOrder.end(), sortVariablesByOrder);

            /**
             * reorder variables
             */
            std::vector<size_t> varIndexOld2New(varInfo.size(), -1);
            for (size_t j = 0; j < varOrder.size(); ++j) {
                varIndexOld2New[varOrder[j].originalIndex] = j;
            }

            newVarInfo.resize(varInfo.size());
            for (size_t j = 0; j < varOrder.size(); ++j) {
                newVarInfo[j] = varInfo[varOrder[j].originalIndex];
                int oldDerivOfIndex = newVarInfo[j].getAntiDerivative();
                if (oldDerivOfIndex >= 0)
                    newVarInfo[j].setAntiDerivative(varIndexOld2New[oldDerivOfIndex]);
                int oldDerivIndex = newVarInfo[j].getDerivative();
                if (oldDerivIndex >= 0)
                    newVarInfo[j].setDerivative(varIndexOld2New[oldDerivIndex]);
            }

            /**
             * reorder equations
             */
            newEqInfo = eqInfo; //copy
            for (size_t i = 0; i < newEqInfo.size(); i++) {
                int oldVIndex = newEqInfo[i].getAssignedVarIndex();
                if (oldVIndex >= 0) {
                    newEqInfo[i].setAssignedVarIndex(varIndexOld2New[oldVIndex]);
                }
            }

            // sort by the order of the assigned variables
            std::sort(newEqInfo.begin(), newEqInfo.end(), sortEquationByAssignedOrder);

            /**
             * Generate an operation graph
             */
            CodeHandler<Base> handler;

            vector<CGBase> indep0(fun->Domain());
            handler.makeVariables(indep0);

            const vector<CGBase> res0 = fun->Forward(0, indep0);

            /**
             * Implement the reordering in the model
             */
            std::auto_ptr<ADFun<CGBase > > reorderedFun(generateReorderedModel(handler, res0, varInfo, newVarInfo, eqInfo, newEqInfo));

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "reordered DAE equations and variables:\n";
            Plantelides<Base>::printModel(reorderedFun.get(), newVarInfo);
#endif

            return reorderedFun;
        }

        inline ADFun<CGBase>* generateReorderedModel(CodeHandler<Base>& handler,
                                                     const std::vector<CGBase>& res0,
                                                     const std::vector<DaeVarInfo>& varInfo,
                                                     const std::vector<DaeVarInfo>& newVarInfo,
                                                     const std::vector<DaeEquationInfo>& eqInfo,
                                                     const std::vector<DaeEquationInfo>& newEqInfo) const {
            using std::vector;

            vector<ADCG> indepNewOrder(handler.getIndependentVariableSize());
            assert(indepNewOrder.size() == newVarInfo.size());

            for (size_t p = 0; p < newVarInfo.size(); p++) {
                int origIndex = newVarInfo[p].getOriginalIndex();
                if (origIndex >= 0) {
                    indepNewOrder[p] = this->x_[origIndex];
                }
            }

            Independent(indepNewOrder);

            /**
             * the model must be called with the handler order
             * 
             * removed variables using substitution are taken out from the list
             * of independent variables in the handler
             */
            std::set<size_t> newVarIndexes;
            for (size_t j = 0; j < newVarInfo.size(); j++) {
                newVarIndexes.insert(newVarInfo[j].getOriginalIndex());
            }

            std::map<size_t, size_t> varOrig2HandlerIndex;
            size_t handlerIndex = 0;
            for (size_t j = 0; j < varInfo.size(); j++) {
                int orig = varInfo[j].getOriginalIndex();
                if (newVarIndexes.find(orig) != newVarIndexes.end()) {
                    varOrig2HandlerIndex[orig] = handlerIndex++;
                }
            }

            vector<ADCG> indepHandlerOrder(handler.getIndependentVariableSize());
            for (size_t p = 0; p < newVarInfo.size(); p++) {
                size_t origIndex = newVarInfo[p].getOriginalIndex();
                indepHandlerOrder[varOrig2HandlerIndex[origIndex]] = indepNewOrder[p];
            }

            // reorder equations
            std::map<size_t, size_t> eqOrigIndex2OldIndex;
            for (size_t i = 0; i < eqInfo.size(); i++) {
                eqOrigIndex2OldIndex[eqInfo[i].getOriginalIndex()] = i;
            }

            vector<CGBase> resNewOrder(newEqInfo.size());
            for (size_t i = 0; i < newEqInfo.size(); i++) {
                size_t oldIndex = eqOrigIndex2OldIndex[newEqInfo[i].getOriginalIndex()];
                resNewOrder[i] = res0[oldIndex];
            }

            // evaluate the model
            Evaluator<Base, CGBase> evaluator0(handler, resNewOrder);
            vector<ADCG> depNewOrder = evaluator0.evaluate(indepHandlerOrder);

            return new ADFun<CGBase > (indepNewOrder, depNewOrder);
        }

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
                for (size_t j = diffVarStart_; j < this->vnodes_.size(); j++) {
                    assert(this->vnodes_[j]->antiDerivative() != NULL);
                    size_t t = this->vnodes_[j]->tapeIndex();
                    if (jacSparsity_[i * n + t]) {
                        row.push_back(i);
                        col.push_back(t);
                    }
                }
            }

            vector<CG<Base> > jac(row.size());

            vector<CG<Base> > indep(n);
            std::copy(this->x_.begin(), this->x_.end(), indep.begin());
            std::fill(indep.begin() + this->x_.size(), indep.end(), 0);

            CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
            this->reducedFun_->SparseJacobianReverse(indep, jacSparsity_,
                                                     row, col, jac, work);

            // resize and zero matrix
            jacobian_.resize(m - diffEqStart_, this->vnodes_.size() - diffVarStart_);

            map<size_t, Vnode<Base>*> origIndex2var;
            for (size_t j = diffVarStart_; j< this->vnodes_.size(); j++) {
                Vnode<Base>* jj = this->vnodes_[j];
                origIndex2var[jj->tapeIndex()] = jj;
            }

            // normalize values
            for (size_t e = 0; e < jac.size(); e++) {
                Enode<Base>* eqOrig = this->enodes_[row[e]]->originalEquation();
                Vnode<Base>* vOrig = origIndex2var[col[e]]->originalVariable(this->origTimeDependentCount_);

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

        /**
         * Determines the best matrix
         * \param mat The matrix
         * @return The best condition value (lowest possible and real)
         */
        inline static Base evalBestMatrixCondition(const MatrixB& mat) {

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

        inline static void printModel(CodeHandler<Base>& handler,
                                      const std::vector<CGBase>& res,
                                      const std::vector<DaeVarInfo>& varInfo,
                                      const std::set<size_t>& erasedVariables,
                                      const std::set<size_t>& erasedEquations) {
            std::vector<std::string> indepNames;
            for (size_t p = 0; p < varInfo.size(); p++) {
                if (erasedVariables.find(p) == erasedVariables.end()) {
                    // not erased from model
                    indepNames.push_back(varInfo[p].getName());
                }
            }
            assert(handler.getIndependentVariableSize() == indepNames.size());

            CLanguage<Base> lang("double");
            std::vector<CGBase> resAux;
            for (size_t p = 0; p < res.size(); ++p) {
                if (erasedEquations.find(p) == erasedEquations.end()) {
                    resAux.push_back(res[p]);
                }
            }
            std::vector<std::string> depNames;
            CLangCustomVariableNameGenerator<Base> nameGen(depNames, indepNames);
            handler.generateCode(std::cout, lang, resAux, nameGen);
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
