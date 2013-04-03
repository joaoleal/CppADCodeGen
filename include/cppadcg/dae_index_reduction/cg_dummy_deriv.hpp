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

            std::auto_ptr<ADFun<CG<Base> > > fun(Plantelides<Base>::reduceIndex(reducedVarInfo, reducedEqInfo));
            if (fun.get() == NULL)
                return NULL; //nothing to do (no index reduction required)

            newEqInfo = reducedEqInfo; // copy
            addDummyDerivatives(reducedVarInfo, newVarInfo);

            if (reduceEquations_) {
                std::vector<DaeVarInfo> varInfo = newVarInfo; // copy
                std::auto_ptr<ADFun<CG<Base> > > funShort = reduceEquations(varInfo, newVarInfo,
                                                                            reducedEqInfo, newEqInfo);
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
                assert((*j)->tapeIndex() >= 0);
                assert((*j)->antiDerivative() != NULL);
                assert((*j)->antiDerivative()->tapeIndex() >= 0);

                newVarInfo[(*j)->tapeIndex()].setAntiDerivative(-1);
                newVarInfo[(*j)->antiDerivative()->tapeIndex()].setDerivative(-1);
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

            assert(reducedVarInfo.size() == this->reducedFun_->Domain());
            assert(reducedEqInfo.size() == this->reducedFun_->Range());

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
#ifdef CPPAD_CG_DAE_VERBOSE
                    std::cout << "unable to use equation " << i->index() << " to solve for variable " << dummy->name() << ": " << ex.what() << std::endl;
#endif
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
            typedef vector<SourceCodePathNode<Base> > SourceCodePath;

            newEqInfo = eqInfo; // copy (we will have the same number of equations)

            /**
             * Generate an operation graph
             */
            CodeHandler<Base> handler;

            vector<CGBase> indep0(fun.Domain());
            handler.makeVariables(indep0);

            vector<CGBase> res0 = fun.Forward(0, indep0);

            vector<bool> jacSparsity = jacobianSparsity < vector<bool> >(fun);

            std::map<int, int> assignedVar2Eq;
            for (size_t i = 0; i < eqInfo.size(); ++i) {
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
                 * (the assigned equation might not be the best)
                 */
                SourceCodePath bestPath;
                int bestEquation = -1;
                for (size_t i = 0; i < eqInfo.size(); ++i) {
                    DaeEquationInfo& newEq = newEqInfo[i];

                    bool exists = jacSparsity[varInfo.size() * i + j];
                    bool forOtherDiffVar;
                    if (newEq.getAssignedVarIndex() < 0) {
                        forOtherDiffVar = false;
                    } else {
                        assert(newEq.getAssignedVarIndex() < int(varInfo.size()));
                        const DaeVarInfo& assigned = varInfo[newEq.getAssignedVarIndex()];
                        forOtherDiffVar = (assigned.getAntiDerivative() != jj.getAntiDerivative() && assigned.getAntiDerivative() >= 0);
                    }

                    if (exists && // the variable is present in this equation
                            !newEq.isExplicit() && // the equation was not used before for a different time derivative
                            !forOtherDiffVar// it is not destined for a different time derivative
                            ) {
                        CGBase& dep = res0[i]; // the equation residual

                        vector<SourceCodePath> paths = findPaths(*dep.getSourceCodeFragment(), *indep.getSourceCodeFragment(), 2);
                        if (paths.size() == 1) {
                            if (bestPath.empty() || bestPath.size() > paths[0].size()) {
                                bestPath = paths[0];
                                bestEquation = i;
                            }
                        }
                    }
                }
                if (bestEquation == -1) {
                    throw CGException("Failed to generate semi-explicit DAE: unable to create an explicit equation for " + jj.getName());
                }
                try {
                    CGBase& dep = res0[bestEquation]; // the equation residual

                    handler.substituteIndependent(indep, dep); // removes indep from the list of variables

                    SourceCodeFragment<Base>* alias = indep.getSourceCodeFragment();
                    assert(alias != NULL && alias->operation() == CGAliasOp);
                    dep.getSourceCodeFragment()->makeAlias(alias->arguments()[0]);

                    // it is now an explicit differential equation
                    newEqInfo[bestEquation].setExplicit(true);
                    // possibly switch assigned variables
                    if (newEqInfo[bestEquation].getAssignedVarIndex() != int(j)) {
                        DaeEquationInfo& assignedEq = newEqInfo[assignedVar2Eq.at(j)];
                        assignedEq.setAssignedVarIndex(newEqInfo[bestEquation].getAssignedVarIndex());
                    }
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

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "Semi-Eplicit DAE:\n";
            Plantelides<Base>::printModel(semiExplicitFun.get(), newVarInfo);
#endif

            return semiExplicitFun;
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
                                           MatrixB& work) throw (CGException) {

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
#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "subset Jac:\n" << work << "\n";
#endif

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

            std::set<Vnode<Base>* > newDummies;
            for (int i = 0; i < work.rows(); i++) {
                newDummies.insert(varsLocal[indices(i)]);
            }

#ifdef CPPAD_CG_DAE_VERBOSE
            std::cout << "## new dummy derivatives: "; //"(condition = " << bestCond << "): ";
            for (typename std::set<Vnode<Base>* >::const_iterator it = newDummies.begin(); it != newDummies.end(); ++it)
                std::cout << **it << "; ";
            std::cout << " \n\n";
#endif

            dummyD_.insert(newDummies.begin(), newDummies.end());
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

    };
}

#endif
