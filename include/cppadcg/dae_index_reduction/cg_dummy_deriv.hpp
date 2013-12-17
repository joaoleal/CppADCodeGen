#ifndef CPPAD_CG_DUMMY_DERIV_INCLUDED
#define CPPAD_CG_DUMMY_DERIV_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/QR>

#include <cppadcg/dae_index_reduction/cg_pantelides.hpp>
//#include <unsupported/Eigen/NonLinearOptimization>

namespace CppAD {

    /**
     * Sorts variable nodes according to the variable differentiation order
     * 
     * @param i
     * @param j
     * @return true if i should come before j
     */
    template<class Base>
    bool sortVnodesByOrder(Vnode<Base>* i, Vnode<Base>* j) {
        return (i->order() > j->order());
    }

    /**
     * Utility class used to sort varibles in the DAE
     */
    class DaeVarOrderInfo {
    public:
        size_t originalIndex;
        size_t originalIndex0;
        bool hasDerivatives;
        int order;

        inline DaeVarOrderInfo() :
            originalIndex(0),
            originalIndex0(0),
            hasDerivatives(false),
            order(-1) {
        }

        inline DaeVarOrderInfo(size_t moriginalIndex, size_t moriginalIndex0, bool mhasDerivatives, int morder) :
            originalIndex(moriginalIndex),
            originalIndex0(moriginalIndex0),
            hasDerivatives(mhasDerivatives),
            order(morder) {
        }
    };

    /**
     * Utility class used to sort equations in the DAE system
     */
    class DaeEqOrderInfo {
    public:
        size_t originalIndex;
        size_t originalIndex0;
        bool differential;
        int assignedVar;

        inline DaeEqOrderInfo() :
            originalIndex(0),
            originalIndex0(0),
            differential(false),
            assignedVar(-1) {
        }

        inline DaeEqOrderInfo(size_t moriginalIndex, size_t moriginalIndex0, bool mdifferential, int massignedVar) :
            originalIndex(moriginalIndex),
            originalIndex0(moriginalIndex0),
            differential(mdifferential),
            assignedVar(massignedVar) {
        }
    };

    /**
     * Sorts variables based on the differentiation order, whether they are 
     * algebraic or differential and the order in the original model
     * 
     * @param i
     * @param j
     * @return true if i should come before j
     */
    inline bool sortVariablesByOrder(const DaeVarOrderInfo& i, const DaeVarOrderInfo& j) {
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
     * Sorts equations according to the equation type (differential/algebraic)
     * and original index
     * 
     * @param i
     * @param j
     * @return true if i should come before j
     */
    inline bool sortEquationByAssignedOrder2(const DaeEqOrderInfo& i, const DaeEqOrderInfo& j) {
        if (i.differential) {
            if (j.differential)
                return i.assignedVar < j.assignedVar;
            else
                return true;
        } else {
            if (j.differential) {
                return false;
            } else {
                if (i.originalIndex0 == j.originalIndex0) {
                    return i.originalIndex == j.originalIndex0;
                } else {
                    return i.originalIndex0 < j.originalIndex0;
                }
            }
        }
    }

    /**
     * Dummy derivatives DAE index reduction algorithm
     */
    template<class Base>
    class DummyDerivatives : public Pantelides<Base> {
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
        std::vector<Vnode<Base>* > dummyD_;
        bool reduceEquations_;
        bool generateSemiExplicitDae_;
        bool reorder_;
    public:

        /**
         * Creates the DAE index reduction algorithm that implements the dummy
         * derivatives method.
         * 
         * @param fun The DAE model
         * @param varInfo DAE model variable classification
         * @param x typical variable values (used to determine Jacobian values)
         * @param normVar variable normalization values
         * @param normEq equation normalization values
         */
        DummyDerivatives(ADFun<CG<Base> >* fun,
                         const std::vector<DaeVarInfo>& varInfo,
                         const std::vector<Base>& x,
                         const std::vector<Base>& normVar,
                         const std::vector<Base>& normEq) :
            Pantelides<Base>(fun, varInfo, x),
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
         * @param newVarInfo Variable related information of the reduced index
         *                   model
         * @param equationInfo Equation related information of the reduced index
         *                     model
         * @return the reduced index model (must be deleted by user)
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

            std::auto_ptr<ADFun<CG<Base> > > fun(Pantelides<Base>::reduceIndex(reducedVarInfo, reducedEqInfo));
            if (fun.get() == NULL)
                return NULL; //nothing to do (no index reduction required)

            newEqInfo = reducedEqInfo; // copy
            addDummyDerivatives(reducedVarInfo, newVarInfo);

            if (reduceEquations_ || generateSemiExplicitDae_) {

                matchVars2Eqs4Eliminaton(newVarInfo, newEqInfo);

                if (reduceEquations_) {
                    std::vector<DaeVarInfo> varInfo = newVarInfo; // copy
                    std::vector<DaeEquationInfo> eqInfo = newEqInfo; // copy
                    std::auto_ptr<ADFun<CG<Base> > > funShort = reduceEquations(varInfo, newVarInfo,
                                                                                eqInfo, newEqInfo);
                    fun = funShort;
                }

                if (generateSemiExplicitDae_) {
                    std::vector<DaeVarInfo> varInfo = newVarInfo; // copy
                    std::vector<DaeEquationInfo> eqInfo = newEqInfo; // copy
                    std::auto_ptr<ADFun<CG<Base> > > semiExplicit = generateSemiExplicitDAE(*fun.get(),
                                                                                            varInfo, newVarInfo,
                                                                                            eqInfo, newEqInfo);
                    fun = semiExplicit;
                }
            }

            if (reorder_) {
                std::vector<DaeVarInfo> varInfo = newVarInfo; // copy
                std::vector<DaeEquationInfo> eqInfo = newEqInfo; // copy
                std::auto_ptr<ADFun<CG<Base> > > reorderedFun = reorderModelEqNVars(*fun.get(),
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


            MatrixB workJac;

            while (true) {

                if (this->verbosity_ >= VERBOSITY_HIGH) {
                    std::cout << "# equation selection: ";
                    for (size_t i = 0; i < eqs.size(); i++)
                        std::cout << *eqs[i] << "; ";
                    std::cout << "\n";

                    std::cout << "# variable selection: ";
                    for (size_t j = 0; j < vars.size(); j++)
                        std::cout << *vars[j] << "; ";
                    std::cout << "\n";
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
            typename std::vector<Vnode<Base>* >::const_iterator j;
            for (j = dummyD_.begin(); j != dummyD_.end(); ++j) {
                CPPADCG_ASSERT_UNKNOWN((*j)->tapeIndex() >= 0);
                CPPADCG_ASSERT_UNKNOWN((*j)->antiDerivative() != NULL);
                CPPADCG_ASSERT_UNKNOWN((*j)->antiDerivative()->tapeIndex() >= 0);

                newVarInfo[(*j)->tapeIndex()].setAntiDerivative(-1);
                newVarInfo[(*j)->antiDerivative()->tapeIndex()].setDerivative(-1);
            }

            if (this->verbosity_ >= VERBOSITY_HIGH) {
                std::cout << "## dummy derivatives:\n";

                for (j = dummyD_.begin(); j != dummyD_.end(); ++j)
                    std::cout << "# " << **j << "   \t" << newVarInfo[(*j)->tapeIndex()].getName() << "\n";
                std::cout << "# \n";
                Pantelides<Base>::printModel(this->reducedFun_, newVarInfo);
            }

        }

        /**
         * Attempts to reduce the number of equations by variable substitution.
         * 
         * @param newVarInfo Variable information of the resulting model
         * @return The new DAE reduced model with (possibly) less equations and
         *         variables
         */
        inline std::auto_ptr<ADFun<CGBase > > reduceEquations(const std::vector<DaeVarInfo>& reducedVarInfo,
                                                              std::vector<DaeVarInfo>& newVarInfo,
                                                              const std::vector<DaeEquationInfo>& reducedEqInfo,
                                                              std::vector<DaeEquationInfo>& newEqInfo) {
            using namespace std;
            using std::vector;
            using std::map;
            using std::map;

            CPPADCG_ASSERT_UNKNOWN(reducedVarInfo.size() == this->reducedFun_->Domain());
            CPPADCG_ASSERT_UNKNOWN(reducedEqInfo.size() == this->reducedFun_->Range());

            newEqInfo = reducedEqInfo; // copy

            /**
             * Generate an operation graph
             */
            CodeHandler<Base> handler;

            vector<CGBase> indep0(this->reducedFun_->Domain());
            handler.makeVariables(indep0);

            vector<CGBase> res0 = this->reducedFun_->Forward(0, indep0);

            map<int, int> assignedVar2Eq;
            for (size_t i = 0; i < newEqInfo.size(); ++i) {
                DaeEquationInfo& newEq = newEqInfo[i];
                if (newEq.getAssignedVarIndex() >= 0)
                    assignedVar2Eq[newEq.getAssignedVarIndex()] = i;
            }

            /**
             * maps the equations indexes of the reduced model to the new 
             * equation indexes in the model with less equations and variables
             * (removed equations have negative indexes)
             */
            vector<int> eqIndexReduced2Short(this->enodes_.size());
            for (size_t i = 0; i < eqIndexReduced2Short.size(); i++) {
                eqIndexReduced2Short[i] = i;
            }

            /**
             * maps the variables indexes in the tape of the reduced model to 
             * the  new tape indexes in the model with less equations and
             * variables (removed variables have negative indexes)
             */
            vector<int> tapeIndexReduced2Short(reducedVarInfo.size());
            for (size_t j = 0; j < tapeIndexReduced2Short.size(); j++) {
                tapeIndexReduced2Short[j] = j;
            }

            /**
             * attempt to eliminate dummy derivatives and the equation they are
             * assigned to
             */
            set<size_t> erasedVariables;
            set<size_t> erasedEquations;

            typename vector<Vnode<Base>* >::const_iterator j;
            for (j = dummyD_.begin(); j != dummyD_.end(); ++j) {
                Vnode<Base>* dummy = *j;

                /**
                 * Determine which equation to use to eliminate the dummy derivative
                 */
                map<int, int>::const_iterator ita = assignedVar2Eq.find(dummy->tapeIndex());
                if (ita == assignedVar2Eq.end()) {
                    if (this->verbosity_ >= VERBOSITY_HIGH)
                        std::cout << "unable to solve for variable " << dummy->name() << "." << std::endl;

                    continue; // unable to solve for a dummy variable: keep the equation and variable
                }
                int bestEquation = ita->second;

                try {
                    // eliminate all references to the dummy variable by substitution
                    handler.substituteIndependent(indep0[dummy->tapeIndex()], res0[bestEquation]);
                    tapeIndexReduced2Short[dummy->tapeIndex()] = -1;
                    eqIndexReduced2Short[bestEquation] = -1;

                    if (this->verbosity_ >= VERBOSITY_HIGH) {
                        std::cout << "######### use equation " << bestEquation << " to solve for variable " << dummy->name() << std::endl;
                        erasedVariables.insert(dummy->tapeIndex());
                        erasedEquations.insert(bestEquation);
                        printModel(handler, res0, reducedVarInfo, erasedVariables, erasedEquations);
                    }

                } catch (const CGException& ex) {
                    // unable to solve for a dummy variable: keep the equation and variable
                    if (this->verbosity_ >= VERBOSITY_HIGH)
                        std::cout << "unable to use equation " << bestEquation << " to solve for variable " << dummy->name() << ": " << ex.what() << std::endl;
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
            CPPADCG_ASSERT_UNKNOWN(tapeIndexReduced2Short.size() == reducedVarInfo.size());

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

            if (this->verbosity_ >= VERBOSITY_HIGH) {
                std::cout << "DAE with less equations and variables:\n";
                Pantelides<Base>::printModel(shortFun.get(), newVarInfo);
            }

            return shortFun;
        }

        /**
         * Attempts to generate a semi-explicit DAE.
         * 
         * @param reorder place all the differential equations and variables
         *                together
         * @param differentialEqs 
         * @return The new semi-explicit DAE model with less variables (without
         *         the time derivative variables)
         */
        inline std::auto_ptr<ADFun<CGBase > > generateSemiExplicitDAE(ADFun<CG<Base> >& fun,
                                                                      const std::vector<DaeVarInfo>& varInfo,
                                                                      std::vector<DaeVarInfo>& newVarInfo,
                                                                      const std::vector<DaeEquationInfo>& eqInfo,
                                                                      std::vector<DaeEquationInfo>& newEqInfo) throw (CGException) {
            using namespace std;
            using std::vector;
            using std::map;

            newEqInfo = eqInfo; // copy (we will have the same number of equations)

            /**
             * Generate an operation graph
             */
            CodeHandler<Base> handler;

            vector<CGBase> indep0(fun.Domain());
            handler.makeVariables(indep0);

            vector<CGBase> res0 = fun.Forward(0, indep0);

            map<int, int> assignedVar2Eq;
            for (size_t i = 0; i < newEqInfo.size(); ++i) {
                DaeEquationInfo& newEq = newEqInfo[i];
                assignedVar2Eq[newEq.getAssignedVarIndex()] = i;
            }

            /**
             * Eliminate time derivatives from equations
             */
            for (size_t j = 0; j < varInfo.size(); ++j) {
                const DaeVarInfo& jj = varInfo[j];
                if (jj.getAntiDerivative() < 0) {
                    continue; // not a time derivative
                }
                CGBase& indep = indep0[j]; // the time derivative
                /**
                 * Determine which equation to keep as differential
                 */
                map<int, int>::const_iterator ita = assignedVar2Eq.find(j);
                if (ita == assignedVar2Eq.end()) {
                    throw CGException("Failed to generate semi-explicit DAE: unable to create an explicit equation for " + jj.getName());
                }

                int bestEquation = ita->second;

                try {
                    CGBase& dep = res0[bestEquation]; // the equation residual

                    handler.substituteIndependent(indep, dep); // removes indep from the list of variables

                    OperationNode<Base>* alias = indep.getOperationNode();
                    CPPADCG_ASSERT_UNKNOWN(alias != NULL && alias->getOperationType() == CGAliasOp);
                    dep.getOperationNode()->makeAlias(alias->getArguments()[0]); // not a residual anymore but equal to alias (explicit equation)

                    // it is now an explicit differential equation
                    newEqInfo[bestEquation].setExplicit(true);
                    // the derivative variable will disappear, associate the equation with the original variable
                    newEqInfo[bestEquation].setAssignedVarIndex(jj.getAntiDerivative());
                } catch (const CGException& ex) {
                    // unable to solve for a dummy variable: keep the equation and variable
                    throw CGException(string("Failed to generate semi-explicit DAE: ") + ex.what());
                }
            }

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

            for (size_t i = 0; i < newEqInfo.size(); ++i) {
                const DaeEquationInfo& ii = newEqInfo[i];
                int j = ii.getAssignedVarIndex();
                if (j >= 0)
                    newEqInfo[i].setAssignedVarIndex(varIndexOld2New[j]);
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

            if (this->verbosity_ >= VERBOSITY_HIGH) {
                std::cout << "Semi-Eplicit DAE:\n";
                Pantelides<Base>::printModel(semiExplicitFun.get(), newVarInfo);
            }

            return semiExplicitFun;
        }

        inline void matchVars2Eqs4Eliminaton(std::vector<DaeVarInfo>& varInfo,
                                             std::vector<DaeEquationInfo>& eqInfo) throw (CGException) {
            using namespace CppAD::extra;
            using std::vector;
            using std::map;
            typedef vector<OperationPathNode<Base> > SourceCodePath;

            CPPADCG_ASSERT_UNKNOWN(eqInfo.size() == this->enodes_.size());
            CPPADCG_ASSERT_UNKNOWN(varInfo.size() == this->reducedFun_->Domain());
            CPPADCG_ASSERT_UNKNOWN(eqInfo.size() == this->reducedFun_->Range());

            CodeHandler<Base> handler;

            vector<CGBase> indep0(this->reducedFun_->Domain());
            handler.makeVariables(indep0);

            vector<CGBase> res0 = this->reducedFun_->Forward(0, indep0);

            vector<bool> jacSparsity = jacobianSparsity < vector<bool> >(*this->reducedFun_);

            std::map<int, int> origAssignedVar2Eq;
            for (size_t i = 0; i < eqInfo.size(); ++i) {
                if (eqInfo[i].getAssignedVarIndex() >= 0) {
                    origAssignedVar2Eq[eqInfo[i].getAssignedVarIndex()] = i;
                }
            }


            vector<Vnode<Base>*> diffVariables;
            vector<Vnode<Base>*> dummyVariables;
            vector<Vnode<Base>*> variables;
            vector<Enode<Base>*> equations(this->enodes_.size(), NULL);
            try {
                /**
                 * Create a new bipartite graph
                 */
                // create variable nodes
                map<Vnode<Base>*, Vnode<Base>*> eliminateOrig2New;
                for (size_t j = 0; j < this->vnodes_.size(); j++) {
                    Vnode<Base>* v = this->vnodes_[j];
                    if (std::find(dummyD_.begin(), dummyD_.end(), v) != dummyD_.end()) {
                        if (reduceEquations_) {
                            dummyVariables.push_back(new Vnode<Base>(j, v->tapeIndex(), v->name()));
                            eliminateOrig2New[v] = dummyVariables.back();
                        }
                    } else if (v->antiDerivative() != NULL) {
                        if (generateSemiExplicitDae_) {
                            diffVariables.push_back(new Vnode<Base>(j, v->tapeIndex(), v->name()));
                            eliminateOrig2New[v] = diffVariables.back();
                        }
                    }
                }

                variables.reserve(diffVariables.size() + dummyVariables.size());
                variables.insert(variables.end(), diffVariables.begin(), diffVariables.end()); // must be added 1st (priority)
                variables.insert(variables.end(), dummyVariables.begin(), dummyVariables.end());

                // create equation nodes
                for (size_t i = 0; i < this->enodes_.size(); i++) {
                    equations[i] = new Enode<Base>(i);
                    const vector<Vnode<Base>*>& origVars = this->enodes_[i]->originalVariables();
                    for (size_t p = 0; p < origVars.size(); p++) {
                        Vnode<Base>* jOrig = origVars[p];

                        typename map<Vnode<Base>*, Vnode<Base>*>::const_iterator it;
                        it = eliminateOrig2New.find(jOrig);
                        if (it != eliminateOrig2New.end() &&
                                jacSparsity[varInfo.size() * i + jOrig->tapeIndex()]) {
                            Vnode<Base>* j = it->second;

                            CGBase& dep = res0[i]; // the equation residual
                            CGBase& indep = indep0[j->tapeIndex()];

                            vector<SourceCodePath> paths = handler.findPaths(*dep.getOperationNode(), *indep.getOperationNode(), 2);
                            if (paths.size() == 1 && isSolvable(paths[0])) {
                                equations[i]->addVariable(j);
                            }
                        }
                    }
                }

                map<size_t, Vnode<Base>*> tape2FreeVariables;
                typename vector<Vnode<Base>*>::const_iterator itj;
                for (itj = variables.begin(); itj != variables.end(); ++itj) {
                    Vnode<Base>* j = *itj;
                    tape2FreeVariables[j->tapeIndex()] = j;
                }


                /**
                 * Match equations to variables
                 */
                while (true) {
                    size_t assigned;
                    do {
                        do {
                            /**
                             * assign variables that can only be solved by a single equation
                             */
                            do {
                                assigned = 0;
                                typename vector<Vnode<Base>*>::const_iterator itj;
                                for (itj = diffVariables.begin(); itj != diffVariables.end(); ++itj) {
                                    Vnode<Base>* j = *itj;
                                    if (!j->isDeleted() && j->equations().size() == 1) {
                                        Enode<Base>& i = *j->equations()[0];
                                        if (i.assigmentVariable() == NULL) {
                                            if (!assignVar2Equation(i, res0, *j, indep0, handler,
                                                                    jacSparsity, tape2FreeVariables,
                                                                    equations, varInfo)) {
                                                throw CGException("Failed to solve equation " + i.name() + " for variable " + j->name());
                                            }
                                            assigned++;
                                        }
                                    }
                                }
                            } while (assigned > 0);

                            /**
                             * assign dummy derivatives that can only be solved by a single 
                             * equation 
                             */
                            assigned = 0;
                            typename vector<Vnode<Base>*>::const_iterator itj;
                            for (itj = dummyVariables.begin(); itj != dummyVariables.end(); ++itj) {
                                Vnode<Base>* j = *itj;
                                if (!j->isDeleted() && j->equations().size() == 1) {
                                    Enode<Base>& i = *j->equations()[0];
                                    if (i.assigmentVariable() == NULL) {
                                        if (assignVar2Equation(i, res0, *j, indep0, handler,
                                                               jacSparsity, tape2FreeVariables,
                                                               equations, varInfo))
                                            assigned++;
                                    }
                                }
                            }
                        } while (assigned > 0);

                        /**
                         * assign equations that can only be used to solve for 
                         * a single variable
                         */
                        assigned = 0;
                        typename vector<Enode<Base>*>::const_iterator iti;
                        for (iti = equations.begin(); iti != equations.end(); ++iti) {
                            Enode<Base>& i = **iti;
                            if (i.assigmentVariable() == NULL && i.variables().size() == 1) {
                                Vnode<Base>* j = i.variables()[0];
                                if (assignVar2Equation(i, res0, *j, indep0, handler,
                                                       jacSparsity, tape2FreeVariables,
                                                       equations, varInfo))
                                    assigned++;
                            }
                        }

                    } while (assigned > 0);

                    /**
                     * All variables have at least two equations that can be used
                     * and all equations have at least two variables
                     * choose a tearing variable/equation
                     */
                    typename vector<Vnode<Base>*>::const_iterator itj;

                    assigned = 0;
                    for (itj = variables.begin(); itj != variables.end(); ++itj) {
                        Vnode<Base>* j = *itj;
                        if (!j->isDeleted()) {
                            typename vector<Enode<Base>*>::const_iterator iti;
                            for (iti = j->equations().begin(); iti != j->equations().end(); ++iti) {
                                Enode<Base>& i = **iti;
                                if (i.assigmentVariable() == NULL) {
                                    if (assignVar2Equation(i, res0, *j, indep0, handler,
                                                           jacSparsity, tape2FreeVariables,
                                                           equations, varInfo)) {
                                        assigned++;
                                        break;
                                    }
                                }
                            }
                            if (assigned > 0)
                                break;
                        }
                    }

                    if (assigned == 0) {
                        break; // done
                    }
                }

                /**
                 * save results
                 */
                for (itj = variables.begin(); itj != variables.end(); ++itj) {
                    Vnode<Base>* j = *itj;
                    if (j->assigmentEquation() != NULL) {
                        int i = j->assigmentEquation()->index();
                        DaeEquationInfo& eq = eqInfo[i];

                        if (eq.getAssignedVarIndex() != int(j->tapeIndex())) {
                            if (eq.getAssignedVarIndex() >= 0) {
                                // switch variables
                                if (origAssignedVar2Eq.find(j->tapeIndex()) != origAssignedVar2Eq.end()) {
                                    DaeEquationInfo& origAssignedEq = eqInfo[origAssignedVar2Eq.at(j->tapeIndex())];
                                    origAssignedEq.setAssignedVarIndex(eq.getAssignedVarIndex());
                                }
                            }

                            eq.setAssignedVarIndex(j->tapeIndex());
                        }
                    }
                }

                // verify results
                if (generateSemiExplicitDae_) {
                    std::string error;
                    for (itj = diffVariables.begin(); itj != diffVariables.end(); ++itj) {
                        Vnode<Base>* j = *itj;
                        if (j->assigmentEquation() == NULL) {
                            // failed!!!
                            if (!error.empty())
                                error += ",";
                            error += " " + j->name();
                        }
                    }
                    if (!error.empty())
                        throw CGException("Failed to generate semi-explicit DAE. Could not solve system for the following variables:" + error);
                }

            } catch (...) {
                deleteVectorValues(diffVariables);
                deleteVectorValues(dummyVariables);
                deleteVectorValues(equations);
                throw;
            }

            if (this->verbosity_ >= VERBOSITY_HIGH) {
                typename vector<Vnode<Base>*>::const_iterator itj;
                for (itj = variables.begin(); itj != variables.end(); ++itj) {
                    Vnode<Base>* j = *itj;
                    if (j->assigmentEquation() != NULL)
                        std::cout << "## Variable " + j->name() << " assigned to equation " << j->assigmentEquation()->name() << "\n";
                }
                std::cout << std::endl;
            }
            deleteVectorValues(diffVariables);
            deleteVectorValues(dummyVariables);
            deleteVectorValues(equations);
        }

        inline bool assignVar2Equation(Enode<Base>& i, std::vector<CGBase>& res0,
                                       Vnode<Base>& j, std::vector<CGBase>& indep0,
                                       CodeHandler<Base>& handler,
                                       std::vector<bool>& jacSparsity,
                                       const std::map<size_t, Vnode<Base>*>& tape2FreeVariables,
                                       std::vector<Enode<Base>*>& equations,
                                       std::vector<DaeVarInfo>& varInfo) throw (CGException) {
            using namespace std;
            using std::vector;
            using std::map;
            typedef vector<OperationPathNode<Base> > SourceCodePath;

            std::vector<bool> localJacSparsity = jacSparsity;
            const size_t n = varInfo.size();

            /**
             * Implement the assigment in the model
             */
            CGBase& dep = res0[i.index()];
            CGBase& indep = indep0[j.tapeIndex()];

            std::string indepName;
            if (indep.getOperationNode()->getName() != NULL) {
                indepName = *indep.getOperationNode()->getName();
            }


            try {
                handler.substituteIndependent(indep, dep, false); // indep not removed from the list of variables yet

                OperationNode<Base>* alias = indep.getOperationNode();
                CPPADCG_ASSERT_UNKNOWN(alias != NULL && alias->getOperationType() == CGAliasOp);

                // it is now an explicit differential equation
                //newEqInfo[bestEquation].setExplicit(true);
                // the derivative variable will disappear, associate the equation with the original variable
                //newEqInfo[bestEquation].setAssignedVarIndex(jj.getAntiDerivative());
            } catch (const CGException& ex) {
                // unable to solve for a dummy variable: keep the equation and variable
                throw CGException(string("Failed to generate to solve equation ") + i.name() + " for variable " + j.name() + ": " + ex.what());
            }


            /**
             * Update the connections to the other equations affected  by the
             * substitution
             */
            vector<size_t> nnzs;
            typename map<size_t, Vnode<Base>*>::const_iterator it;
            for (it = tape2FreeVariables.begin(); it != tape2FreeVariables.end(); ++it) {
                size_t tapeJ = it->first;
                if (localJacSparsity[n * i.index() + tapeJ] && tapeJ != j.tapeIndex()) {
                    nnzs.push_back(tapeJ);
                }
            }
            set<Enode<Base>*> affected;
            for (size_t e = 0; e < equations.size(); ++e) {
                if (equations[e] != &i && localJacSparsity[n * e + j.tapeIndex()]) {
                    localJacSparsity[n * e + j.tapeIndex()] = false; // eliminated by substitution
                    affected.insert(equations[e]);
                    for (size_t p = 0; p < nnzs.size(); ++p) {
                        localJacSparsity[n * e + nnzs[p]] = true;
                    }
                }
            }

            if (affected.size() > 0) {
                /**
                 * Redetermine solvability for the affected equations
                 */
                map<size_t, set<Enode<Base>*> > solvable;
                for (size_t e = 0; e < equations.size(); ++e) {
                    Enode<Base>* eq = equations[e];
                    if (&i != eq && affected.find(eq) == affected.end()) {
                        // no change
                        for (size_t v = 0; v < eq->variables().size(); ++v) {
                            solvable[eq->variables()[v]->tapeIndex()].insert(eq);
                        }
                    }
                }

                // redetermine solvability
                typename set<Enode<Base>*>::const_iterator itAff;
                for (itAff = affected.begin(); itAff != affected.end(); ++itAff) {
                    Enode<Base>& a = **itAff;
                    for (it = tape2FreeVariables.begin(); it != tape2FreeVariables.end(); ++it) {
                        size_t j = it->first;
                        if (localJacSparsity[n * a.index() + j]) {
                            vector<SourceCodePath> paths = handler.findPaths(*res0[a.index()].getOperationNode(), *indep0[j].getOperationNode(), 2);
                            if (paths.size() == 1 && isSolvable(paths[0])) {
                                solvable[j].insert(&a);
                            }
                        }
                    }
                }

                // check if any variable stops being solvable
                bool ok = true;
                for (it = tape2FreeVariables.begin(); it != tape2FreeVariables.end(); ++it) {
                    Vnode<Base>* v = it->second;
                    if (v == &j)
                        continue;

                    if (v->assigmentEquation() != NULL) {
                        if (affected.count(v->assigmentEquation()) > 0 &&
                                solvable[v->tapeIndex()].count(v->assigmentEquation()) == 0) {
                            ok = false;
                            break;
                        }
                    } else if (solvable[v->tapeIndex()].size() == 0) {
                        ok = false;
                        break;
                    }
                }

                if (!ok) {
                    handler.undoSubstituteIndependent(*indep.getOperationNode());
                    if (indepName.size() > 0) {
                        indep.getOperationNode()->setName(indepName);
                    }
                    return false;
                }

                /**
                 * Implement changes in graph
                 */
                for (itAff = affected.begin(); itAff != affected.end(); ++itAff) {
                    Enode<Base>& a = **itAff;
                    for (it = tape2FreeVariables.begin(); it != tape2FreeVariables.end(); ++it) {
                        size_t v = it->first;
                        if (localJacSparsity[n * a.index() + v]) {
                            if (solvable[v].count(&a) > 0) {
                                a.addVariable(it->second);
                            } else {
                                // not solvable anymore
                                a.deleteNode(it->second);
                            }
                        }
                    }
                }

            }

            /**
             * Finalize model
             */
            handler.removeIndependent(*indep.getOperationNode());

            /**
             * Implement the assignment in the graph
             */
            j.setAssigmentEquation(i, this->verbosity_);
            j.deleteNode(this->verbosity_);

            jacSparsity = localJacSparsity;

            return true;
        }

        inline std::auto_ptr<ADFun<CGBase > > reorderModelEqNVars(ADFun<CG<Base> >& fun,
                                                                  const std::vector<DaeVarInfo>& varInfo,
                                                                  std::vector<DaeVarInfo>& newVarInfo,
                                                                  const std::vector<DaeEquationInfo>& eqInfo,
                                                                  std::vector<DaeEquationInfo>& newEqInfo) {

            using namespace std;
            using std::vector;

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
            std::vector<DaeVarOrderInfo> varOrder(varInfo.size());
            for (size_t j = 0; j < varInfo.size(); j++) {
                size_t j0;
                int derivOrder = this->determineVariableDiffOrder(varInfo, j, j0);
                if (varInfo[j].isIntegratedVariable()) {
                    derivOrder = -2; // so that it goes last
                }
                bool hasDerivatives = oldVarWithDerivatives.find(j) != oldVarWithDerivatives.end();
                varOrder[j] = DaeVarOrderInfo(j, j0, hasDerivatives, derivOrder);
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

            std::vector<DaeEqOrderInfo> eqOrder(newEqInfo.size());
            for (size_t i = 0; i < newEqInfo.size(); i++) {
                int assignedVar = newEqInfo[i].getAssignedVarIndex();
                size_t i0 = i;
                while (newEqInfo[i0].getAntiDerivative() >= 0) {
                    i0 = newEqInfo[i0].getAntiDerivative();
                }
                bool isDifferential = newEqInfo[i].isExplicit() || (assignedVar >= 0 && newVarInfo[assignedVar].getAntiDerivative() >= 0);
                eqOrder[i] = DaeEqOrderInfo(i, i0, isDifferential, assignedVar);
            }

            std::sort(eqOrder.begin(), eqOrder.end(), sortEquationByAssignedOrder2);

            std::vector<DaeEquationInfo> newEqInfo2(newEqInfo.size());
            for (size_t i = 0; i < eqOrder.size(); i++) {
                newEqInfo2[i] = newEqInfo[eqOrder[i].originalIndex];
            }
            newEqInfo = newEqInfo2;


            /**
             * Generate an operation graph
             */
            CodeHandler<Base> handler;

            vector<CGBase> indep0(fun.Domain());
            handler.makeVariables(indep0);

            const vector<CGBase> res0 = fun.Forward(0, indep0);

            /**
             * Implement the reordering in the model
             */
            std::auto_ptr<ADFun<CGBase > > reorderedFun(generateReorderedModel(handler, res0, varInfo, newVarInfo, eqInfo, newEqInfo));

            if (this->verbosity_ >= VERBOSITY_HIGH) {
                std::cout << "reordered DAE equations and variables:\n";
                Pantelides<Base>::printModel(reorderedFun.get(), newVarInfo);
            }

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
            CPPADCG_ASSERT_UNKNOWN(indepNewOrder.size() == newVarInfo.size());

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
            std::set<size_t> newIds;
            for (size_t j = 0; j < newVarInfo.size(); j++) {
                newIds.insert(newVarInfo[j].getId());
            }

            std::map<size_t, size_t> varId2HandlerIndex;
            size_t handlerIndex = 0; // start the variable count again since some variable might have been removed
            for (size_t j = 0; j < varInfo.size(); j++) {
                int id = varInfo[j].getId();
                if (newIds.find(id) != newIds.end()) {
                    varId2HandlerIndex[id] = handlerIndex++; // not removed from model
                }
            }

            vector<ADCG> indepHandlerOrder(handler.getIndependentVariableSize());
            for (size_t p = 0; p < newVarInfo.size(); p++) {
                size_t id = newVarInfo[p].getId();
                indepHandlerOrder[varId2HandlerIndex[id]] = indepNewOrder[p];
            }

            // reorder equations
            std::map<size_t, size_t> eqId2OldIndex;
            for (size_t i = 0; i < eqInfo.size(); i++) {
                eqId2OldIndex[eqInfo[i].getId()] = i;
            }

            vector<CGBase> resNewOrder(newEqInfo.size());
            for (size_t i = 0; i < newEqInfo.size(); i++) {
                size_t oldIndex = eqId2OldIndex[newEqInfo[i].getId()];
                resNewOrder[i] = res0[oldIndex];
            }

            // evaluate the model
            Evaluator<Base, CGBase> evaluator0(handler, resNewOrder);
            vector<ADCG> depNewOrder = evaluator0.evaluate(indepHandlerOrder);

            return new ADFun<CGBase > (indepNewOrder, depNewOrder);
        }

        /**
         * Determines the Jacobian relative to the differential variables
         * (e.g. dxdt)
         */
        inline void determineJacobian() {
            using namespace CppAD::extra;
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
                    CPPADCG_ASSERT_UNKNOWN(this->vnodes_[j]->antiDerivative() != NULL);
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
                Base normVal = jac[e].getValue() * normVar_[vOrig->tapeIndex()]
                        / normEq_[eqOrig->index()];

                size_t i = row[e]; // same order
                size_t j = origIndex2var[col[e]]->index(); // different order than in model/tape

                jacobian_.coeffRef(i - diffEqStart_, j - diffVarStart_) = normVal;
            }

            jacobian_.makeCompressed();

            if (this->verbosity_ >= VERBOSITY_HIGH) {
                cout << "partial jacobian:\n" << jacobian_ << "\n\n";
                //cout << jacobian_.triangularView<Eigen::Lower > () << "\n\n";
            }
        }

        inline void selectDummyDerivatives(const std::vector<Enode<Base>* >& eqs,
                                           const std::vector<Vnode<Base>* >& vars,
                                           MatrixB& work) throw (CGException) {

            if (eqs.size() == vars.size()) {
                dummyD_.insert(dummyD_.end(), vars.begin(), vars.end());
                if (this->verbosity_ >= VERBOSITY_HIGH) {
                    std::cout << "# new dummy derivatives: ";
                    for (size_t j = 0; j < vars.size(); j++)
                        std::cout << *vars[j] << "; ";
                    std::cout << " \n";
                }
#ifndef NDEBUG
                for (typename std::vector<Vnode<Base>* >::const_iterator it = vars.begin(); it != vars.end(); ++it) {
                    CPPADCG_ASSERT_UNKNOWN(std::find(dummyD_.begin(), dummyD_.end(), *it) == dummyD_.end());
                }
#endif
                return;
            }

            /**
             * Determine the columns that must be removed
             */
            std::set<size_t> excludeCols;
            for (size_t j = 0; j < vars.size(); j++) {
                Vnode<Base>* jj = vars[j];
                bool notZero = false;
                for (size_t i = 0; i < eqs.size(); i++) {
                    Enode<Base>* ii = eqs[i];
                    Base val = jacobian_.coeff(ii->index() - diffEqStart_, jj->index() - diffVarStart_);
                    if (val != Base(0.0)) {
                        notZero = true;
                        break;
                    }
                }
                if (!notZero) {
                    // all zeros: must not choose this column/variable
                    excludeCols.insert(j);
                }
            }

            std::vector<Vnode<Base>* > varsLocal;
            varsLocal.reserve(vars.size() - excludeCols.size());
            for (size_t j = 0; j < vars.size(); j++) {
                if (excludeCols.find(j) == excludeCols.end()) {
                    varsLocal.push_back(vars[j]);
                }
            }


            work.setZero(eqs.size(), varsLocal.size());

            // determine the rows that only contain a single nonzero (a single column)
            for (size_t i = 0; i < eqs.size(); i++) {
                Enode<Base>* ii = eqs[i];
                for (size_t j = 0; j < varsLocal.size(); j++) {
                    Vnode<Base>* jj = varsLocal[j];
                    Base val = jacobian_.coeff(ii->index() - diffEqStart_, jj->index() - diffVarStart_);
                    if (val != Base(0.0)) {
                        work(i, j) = val;
                    }
                }
            }

            if (this->verbosity_ >= VERBOSITY_HIGH)
                std::cout << "subset Jac:\n" << work << "\n";

            Eigen::ColPivHouseholderQR<MatrixB> qr(work);
            qr.compute(work);
            if (qr.rank() < work.rows()) {
                throw CGException("Failed to select dummy derivatives! "
                                  "The resulting system is probably singular for the provided data.");
            }

            typedef typename Eigen::ColPivHouseholderQR<MatrixB>::PermutationType PermutationMatrix;
            typedef typename PermutationMatrix::IndicesType Indices;

            const PermutationMatrix& p = qr.colsPermutation();
            const Indices& indices = p.indices();
            if (indices.size() < work.rows()) {
                throw CGException("Failed to select dummy derivatives! "
                                  "The resulting system is probably singular for the provided data.");
            }

            std::vector<Vnode<Base>* > newDummies;
            for (int i = 0; i < work.rows(); i++) {
                newDummies.push_back(varsLocal[indices(i)]);
            }

            if (this->verbosity_ >= VERBOSITY_HIGH) {
                std::cout << "## new dummy derivatives: "; //"(condition = " << bestCond << "): ";
                for (typename std::vector<Vnode<Base>* >::const_iterator it = newDummies.begin(); it != newDummies.end(); ++it)
                    std::cout << **it << "; ";
                std::cout << " \n\n";
            }
#ifndef NDEBUG
            for (typename std::vector<Vnode<Base>* >::const_iterator it = newDummies.begin(); it != newDummies.end(); ++it) {
                CPPADCG_ASSERT_UNKNOWN(std::find(dummyD_.begin(), dummyD_.end(), *it) == dummyD_.end());
            }
#endif

            dummyD_.insert(dummyD_.end(), newDummies.begin(), newDummies.end());
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
            CPPADCG_ASSERT_UNKNOWN(handler.getIndependentVariableSize() == indepNames.size());

            CLanguage<Base> lang("double");
            vector<CGBase> resAux;
            for (size_t p = 0; p < res.size(); ++p) {
                if (erasedEquations.find(p) == erasedEquations.end()) {
                    resAux.push_back(res[p]);
                }
            }
            std::vector<std::string> depNames;
            CLangCustomVariableNameGenerator<Base> nameGen(depNames, indepNames);
            handler.generateCode(std::cout, lang, resAux, nameGen);
        }

        inline static void printGraphSparsity(const std::vector<bool>& jacSparsity,
                                              const std::map<size_t, Vnode<Base>*>& tape2FreeVariables,
                                              const std::vector<Enode<Base>*>& equations,
                                              const size_t n) {
            for (size_t e = 0; e < equations.size(); ++e) {
                Enode<Base>* eq = equations[e];
                size_t count = 0;
                typename std::map<size_t, Vnode<Base>*>::const_iterator it;
                for (it = tape2FreeVariables.begin(); it != tape2FreeVariables.end(); ++it) {
                    if (jacSparsity[n * eq->index() + it->first]) {
                        if (count == 0)
                            std::cout << "# Equation " << e << ": \t";
                        std::cout << " " << it->second->name();
                        count++;
                    }
                }
                if (count > 0)
                    std::cout << "\n";
            }

            std::cout << std::endl;
        }

        template<class T>
        inline static void deleteVectorValues(std::vector<T*>& v) {
            for (size_t i = 0; i < v.size(); i++) {
                delete v[i];
            }
            v.clear();
        }

    };
}

#endif
