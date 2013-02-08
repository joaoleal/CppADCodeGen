#ifndef CPPAD_CG_PANTELIDES_INCLUDED
#define CPPAD_CG_PANTELIDES_INCLUDED
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

#include <cppadcg/dae_index_reduction/cg_bipartite.hpp>
#include <cppadcg/dae_index_reduction/cg_dae_index_reduction.hpp>
#include <cppadcg/dae_index_reduction/cg_dae_equation_info.hpp>
#include <cppadcg/dae_index_reduction/cg_time_diff.hpp>

namespace CppAD {

    /**
     * Pantelides DAE index reduction algorithm
     */
    template<class Base>
    class Plantelides : public DaeIndexReduction<Base> {
    protected:
        typedef CG<Base> CGBase;
        typedef AD<CGBase> ADCG;
    protected:
        // typical values;
        std::vector<Base> x_;
        // original sparsity
        std::vector<bool> sparsity_;
        // Bipartite graph ([equation i][variable j])
        std::vector<Vnode<Base>*> vnodes_;
        std::vector<Enode<Base>*> enodes_;
        // new index reduced model
        ADFun<CGBase>* reducedFun_;
        // the maximum order of the time derivatives in the original model
        int origMaxTimeDivOrder_;
    private:
        size_t timeOrigVarIndex_; // time index in the original user model (may not exist)
        size_t timeVarIndex_; // the time index in the graph (may not exist)
    public:

        /**
         * Creates the DAE index reduction algorithm that implements the 
         * Pantelides method.
         * 
         * \param fun The DAE model
         * \param varInfo DAE model variable classification
         * \param x typical variable values (used to avoid NaNs in CppAD checks)
         */
        Plantelides(ADFun<CG<Base> >* fun,
                    const std::vector<DaeVarInfo>& varInfo,
                    const std::vector<Base>& x) :
            DaeIndexReduction<Base>(fun, varInfo),
            x_(x),
            reducedFun_(NULL),
            origMaxTimeDivOrder_(0),
            timeOrigVarIndex_(fun->Domain()),
            timeVarIndex_(fun->Domain()) {

            using namespace std;
            using std::vector;

            const size_t m = fun->Range(); // equation count
            const size_t n = fun->Domain(); // total variable count

            // create equation nodes
            enodes_.reserve(1.2 * m + 1);
            enodes_.resize(m);
            for (size_t i = 0; i < m; i++) {
                enodes_[i] = new Enode<Base > (i);
            }

            // create the variable nodes
            vnodes_.reserve(1.2 * n + 1);
            vnodes_.resize(n);

            // locate the time variable (if present)
            for (size_t dj = 0; dj < n; dj++) {
                if (varInfo[dj].isIntegratedVariable()) {
                    if (timeOrigVarIndex_ < n) {
                        throw CGException("More than one time variable (integrated variable) defined");
                    }
                    timeOrigVarIndex_ = dj;
                }
            }

            // determine the order of each time derivative
            vector<int> derivOrder = determineVariableDiffOrder(varInfo);

            origMaxTimeDivOrder_ = *std::max_element(derivOrder.begin(), derivOrder.end());

            // sort the variables according to the time derivative order
            vector<size_t> new2Orig;
            vector<size_t> orig2New(n);
            new2Orig.reserve(n);
            for (int order = -1; order <= origMaxTimeDivOrder_; order++) {
                for (size_t j = 0; j < n; j++) {
                    if (derivOrder[j] == order) {
                        orig2New[j] = new2Orig.size();
                        new2Orig.push_back(j);
                    }
                }
            }

            std::string timeVarName;
            if (timeOrigVarIndex_ == fun->Domain() || varInfo[timeOrigVarIndex_].getName().empty()) {
                timeVarName = "t";
            } else {
                timeVarName = varInfo[timeOrigVarIndex_].getName();
            }

            stringstream ss;
            size_t paramCount = 0;
            size_t timeDepVarCount = 0;
            for (size_t j = 0; j < n; j++) {
                size_t p = new2Orig[j];
                int origIndex = varInfo[p].getDerivativeOf();
                const std::string& customName = varInfo[p].getName();
                if (origIndex < 0) {
                    // generate the variable name
                    if (varInfo[p].isFunctionOfIntegrated()) {
                        if (customName.empty())
                            ss << "x" << timeDepVarCount;
                        else
                            ss << varInfo[p].getName();
                        timeDepVarCount++;
                    } else if (varInfo[p].isIntegratedVariable()) {
                        if (customName.empty())
                            ss << "t";
                        else
                            ss << varInfo[p].getName();
                        timeVarIndex_ = j;
                    } else {
                        if (customName.empty())
                            ss << "p" << paramCount;
                        else
                            ss << varInfo[p].getName();
                        paramCount++;
                    }

                    vnodes_[j] = new Vnode<Base > (j, p, ss.str());
                    ss.str("");
                    ss.clear();

                    if (!varInfo[p].isFunctionOfIntegrated())
                        vnodes_[j]->makeParameter(); // does not depend on time
                } else {
                    Vnode<Base>* derivativeOf = vnodes_[orig2New[origIndex]];
                    std::string name;
                    if (!customName.empty()) {
                        name = customName;
                    } else {
                        name = "d" + derivativeOf->name() + "d" + timeVarName;

                    }
                    vnodes_[j] = new Vnode<Base > (j, p, derivativeOf, name);
                }
            }

            // create the edges
            sparsity_ = jacobianSparsity < vector<bool>, CGBase > (*fun);

            for (size_t i = 0; i < m; i++) {
                for (size_t p = 0; p < n; p++) {
                    size_t j = orig2New[p];
                    if (sparsity_[i * n + p] && !vnodes_[j]->isDeleted()) {
                        enodes_[i]->addVariable(vnodes_[j]);
                    }
                }
            }

            // make sure the system is not under or over determined
            size_t nvar = 0;
            for (size_t j = 0; j < n; j++) {
                const Vnode<Base>* jj = vnodes_[j];
                if (!jj->isParameter() && // exclude constants
                        (jj->derivativeOf() != NULL || // derivatives
                        jj->derivative() == NULL) // algebraic variables
                        ) {
                    nvar++;
                }
            }

            if (nvar != m) {
                stringstream ss;
                ss << "The system is not well determined. "
                        "The of number of equations (" << enodes_.size() << ") does not match the number of unknown variables (" << nvar << ").";
                throw CGException(ss.str());
            }
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
                                                     std::vector<DaeEquationInfo>& equationInfo) throw (CGException) {
            if (reducedFun_ != NULL)
                throw CGException("reduceIndex() can only be called once!");

            detectSubset2Dif();

#ifdef CPPAD_CG_DAE_VERBOSE
            printResultInfo();
#endif

            generateNewModel(newVarInfo, equationInfo);

            return reducedFun_;
        }

        inline void printResultInfo() {
            std::cout << "\nPantelides DAE differentiation index reduction:\n\n"
                    "   Equations count: " << enodes_.size() << "\n";
            typename std::vector<Enode<Base>*>::const_iterator i;
            for (i = enodes_.begin(); i != enodes_.end(); ++i) {
                const Enode<Base>& ii = **i;
                std::cout << "      " << ii.index() << " - " << ii << "\n";
            }

            size_t paramCount = 0;
            typename std::vector<Vnode<Base>*>::const_iterator j;
            for (j = vnodes_.begin(); j != vnodes_.end(); ++j) {
                if ((*j)->isParameter())
                    paramCount++;
            }

            size_t timeIndepCount = vnodes_.size() - paramCount;

            std::cout << "\n   Parameter count: " << paramCount << "\n";
            std::cout << "\n   Variable count: " << timeIndepCount << "\n";
            for (j = vnodes_.begin(); j != vnodes_.end(); ++j) {
                const Vnode<Base>& jj = **j;
                std::cout << "      " << jj.index() << " - " << jj;
                if (jj.assigmentEquation() != NULL) {
                    std::cout << " assigned to " << *jj.assigmentEquation() << "\n";
                } else if (jj.isParameter()) {
                    std::cout << " is a parameter (time independent)\n";
                } else {
                    std::cout << " NOT assigned to any equation\n";
                }
            }

            std::cout << "\n   Degrees of freedom: " << timeIndepCount - enodes_.size() << "\n";
        }

        virtual ~Plantelides() {
            for (size_t i = 0; i < enodes_.size(); i++)
                delete enodes_[i];

            for (size_t j = 0; j < vnodes_.size(); j++)
                delete vnodes_[j];
        }

    protected:

        static inline std::vector<int> determineVariableDiffOrder(const std::vector<DaeVarInfo>& varInfo) {
            size_t n = varInfo.size();
            // determine the order of each time derivative
            std::vector<int> derivOrder(n, 0);
            for (size_t dj = 0; dj < n; dj++) {
                size_t j0;
                derivOrder[dj] = determineVariableDiffOrder(varInfo, dj, j0);
            }

            return derivOrder;
        }

        static inline int determineVariableDiffOrder(const std::vector<DaeVarInfo>& varInfo, size_t index, size_t& j0) {
            int derivOrder = -1;
            j0 = index;
            if (varInfo[index].isFunctionOfIntegrated()) {
                derivOrder = 0;
                while (varInfo[j0].getDerivativeOf() >= 0) {
                    assert(j0 < varInfo.size());
                    assert(varInfo[j0].isFunctionOfIntegrated());
                    derivOrder++;
                    j0 = varInfo[j0].getDerivativeOf();
                }
            }

            return derivOrder;
        }

        /**
         * 
         */
        inline void detectSubset2Dif() {
            Vnode<Base>* jj;
            Enode<Base>* ll;

            size_t Ndash = enodes_.size();
            for (size_t k = 0; k < Ndash; k++) {
                Enode<Base>* i = enodes_[k];
#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << "Outer loop: equation k = " << *i << "\n";
#endif
                bool pathfound = false;
                while (!pathfound) {

                    /**
                     * delete all V-nodes with A!=0 and their incident edges
                     * from the graph
                     */
                    typename std::vector<Vnode<Base>*>::const_iterator j;
                    for (j = vnodes_.begin(); j != vnodes_.end(); ++j) {
                        jj = *j;
                        if (!jj->isDeleted() && jj->derivative() != NULL) {
                            jj->deleteNode();
                        }
                    }

                    uncolorAll();

                    pathfound = augmentPath(*i);

                    if (!pathfound) {
                        const size_t vsize = vnodes_.size(); // the size might change
                        for (size_t l = 0; l < vsize; ++l) {
                            jj = vnodes_[l];
                            if (jj->isColored() && !jj->isDeleted()) {
                                // add new variable derivatives of colored variables
                                Vnode<Base>* jDiff = new Vnode<Base > (vnodes_.size(), vnodes_.size(), jj);
                                vnodes_.push_back(jDiff);
#ifdef CPPAD_CG_DAE_VERBOSE
                                std::cout << "Created " << *jDiff << "\n";
#endif
                            }
                        }

                        const size_t esize = enodes_.size(); // the size might change
                        for (size_t l = 0; l < esize; l++) {
                            ll = enodes_[l];
                            if (ll->isColored()) {
                                // add new derivative equations for colored equations
                                Enode<Base>* lDiff = new Enode<Base > (enodes_.size(), ll);
                                enodes_.push_back(lDiff);

                                // differentiate newI and create edges!!!
                                dirtyDifferentiateEq(*ll, *lDiff);
#ifdef CPPAD_CG_DAE_VERBOSE
                                std::cout << "Created " << *lDiff << "\n";
#endif
                            }
                        }

                        // structural check to avoid infinite recursion
                        for (size_t l = esize; l < enodes_.size(); l++) {
                            ll = enodes_[l];
                            const std::set<Vnode<Base>*>& nvars = ll->originalVariables();
                            bool ok = false;
                            for (typename std::set<Vnode<Base>*>::const_iterator js = nvars.begin(); js != nvars.end(); ++js) {
                                if ((*js)->equations().size() > 1) {
                                    ok = true;
                                    break;
                                }
                            }
                            if (!ok)
                                throw CGException("Invalid equation structure. The model appears to be over-defined.");
                        }

                        for (j = vnodes_.begin(); j != vnodes_.end(); ++j) {
                            jj = *j;
                            if (jj->isColored() && !jj->isDeleted()) {
                                Vnode<Base>* jDiff = jj->derivative();
                                jDiff->setAssigmentEquation(*jj->assigmentEquation()->derivative());
                            }
                        }

                        i = i->derivative();
#ifdef CPPAD_CG_DAE_VERBOSE
                        std::cout << "Set current equation to (i=" << i->index() << ") " << *i << "\n";
#endif
                    }
                }

            }

        }

        /**
         * 
         * \param i The equation node
         * \return true if an augmented path was found
         */
        bool augmentPath(Enode<Base>& i) {
            i.color();

            const std::set<Vnode<Base>*>& vars = i.variables();
            typename std::set<Vnode<Base>*>::const_iterator j;

            // first look for derivative variables
            for (j = vars.begin(); j != vars.end(); ++j) {
                Vnode<Base>* jj = *j;
                if (jj->derivativeOf() != NULL && jj->assigmentEquation() == NULL) {
                    jj->setAssigmentEquation(i);
                    return true;
                }
            }

            // look for algebraic variables
            for (j = vars.begin(); j != vars.end(); ++j) {
                Vnode<Base>* jj = *j;
                if (jj->derivativeOf() == NULL && jj->assigmentEquation() == NULL) {
                    jj->setAssigmentEquation(i);
                    return true;
                }
            }


            for (j = vars.begin(); j != vars.end(); ++j) {
                Vnode<Base>* jj = *j;
                if (!jj->isColored()) {
                    jj->color();
                    Enode<Base>& k = *jj->assigmentEquation();
                    bool pathFound = augmentPath(k);
                    if (pathFound) {
                        jj->setAssigmentEquation(i);
                        return true;
                    }
                }
            }

            return false;
        }

        inline void uncolorAll() {
            typename std::vector<Vnode<Base>*>::const_iterator j;
            for (j = vnodes_.begin(); j != vnodes_.end(); ++j) {
                (*j)->uncolor();
            }
            typename std::vector<Enode<Base>*>::const_iterator i;
            for (i = enodes_.begin(); i != enodes_.end(); ++i) {
                (*i)->uncolor();
            }
        }

        /**
         * Adds a new equation assuming the new equation differential contains 
         * all variables present in the original equation and their time
         * derivatives (not exactly correct but it works because the 
         * potentially extra variables are removed later)
         * 
         * \param i equation node to differentiate
         */
        inline void dirtyDifferentiateEq(Enode<Base>& i, Enode<Base>& newI) throw (CGException) {
            const std::set<Vnode<Base>*>& vars = i.originalVariables();
            typename std::set<Vnode<Base>*>::const_iterator j;
            for (j = vars.begin(); j != vars.end(); ++j) {
                Vnode<Base>* jj = *j;
                newI.addVariable(jj);
                if (jj->derivative() != NULL) {
                    newI.addVariable(jj->derivative());
                }
            }
        }

        /**
         * Creates a new tape for the index 1 model
         */
        inline void generateNewModel(std::vector<DaeVarInfo>& newVarInfo,
                                     std::vector<DaeEquationInfo>& equationInfo) {
            using std::vector;

            vector<vector<Enode<Base>*> > newEquations;

            // find new equations that must be generated by differentiation
            vector<Enode<Base>* > newEqs;
            size_t origM = this->fun_->Range();
            for (size_t i = 0; i < origM; i++) {
                if (enodes_[i]->derivative() != NULL) {
                    assert(enodes_[i]->derivativeOf() == NULL);
                    newEqs.push_back(enodes_[i]->derivative());
                }
            }

            while (newEqs.size() > 0) {
                newEquations.push_back(newEqs);
                newEqs.clear();
                vector<Enode<Base>*>& eqs = newEquations.back();
                for (size_t i = 0; i < eqs.size(); i++) {
                    if (eqs[i]->derivative() != NULL) {
                        newEqs.push_back(eqs[i]->derivative());
                    }
                }
            }

            if (newEquations.empty()) {
                // nothing to do
                return;
            }

            /**
             * Add the relationship between variables and derivatives
             */
            assert(reducedFun_ == NULL);

            {
                CodeHandler<Base> handler;

                vector<CGBase> indep0(this->fun_->Domain());
                handler.makeVariables(indep0);

                const vector<CGBase> dep0 = this->fun_->Forward(0, indep0);

                /**
                 * generate a new tape
                 */
                vector<ADCG> indepNew;
                if (timeOrigVarIndex_ < vnodes_.size()) {
                    indepNew = vector<ADCG > (vnodes_.size()); // variables + time    
                } else {
                    indepNew = vector<ADCG > (vnodes_.size() + 1); // variables + time
                }
                for (size_t j = 0; j < x_.size(); j++) {
                    indepNew[j] = x_[j];
                }
                Independent(indepNew);

                // variables with the relationship between x dxdt and t
                vector<ADCG> indep2 = prepareTimeDependentVariables(indepNew);
                indep2.resize(indep0.size());

                Evaluator<Base, CGBase> evaluator0(handler, dep0);
                vector<ADCG> depNew = evaluator0.evaluate(indep2);
                depNew.resize(enodes_.size());

                try {
                    reducedFun_ = new ADFun<CGBase > (indepNew, depNew);
                } catch (const std::exception& ex) {
                    throw CGException(std::string("Failed to create ADFun: ") + ex.what());
                }


#ifdef CPPAD_CG_DAE_VERBOSE
                printModel(reducedFun_);
#endif
            }


            /**
             * generate the system of equations by repeatedly differentiating
             * and adding equations to the DAE system
             */
            for (size_t d = 0; d < newEquations.size(); d++) {
                vector<Enode<Base>*>& equations = newEquations[d];

                size_t m = reducedFun_->Domain(); // total variable count
                //size_t n = reducedFun_->Range(); // equation count

                /**
                 * register operations from the other equations
                 */
                CodeHandler<Base> handler0;

                vector<CGBase> indep0(m);
                handler0.makeVariables(indep0);

                vector<CGBase> dep = reducedFun_->Forward(0, indep0);

                /**
                 * register operations used to differentiate the equations
                 */
                //forwardTimeDiff(equations, dep);
                reverseTimeDiff(equations, dep);

                delete reducedFun_; // not needed anymore
                reducedFun_ = NULL;

                /**
                 * reconstruct the new system of equations 
                 */
                vector<ADCG> indep2;
                vector<ADCG> indepNew;

                if (d < newEquations.size() - 1) {
                    indepNew.resize(m);
                } else if (timeOrigVarIndex_ == vnodes_.size()) {
                    // the very last model creation
                    indepNew.resize(m - 1); // take out time (it was added by this function and not the user)
                } else {
                    // the very last model creation
                    indepNew.resize(m);
                }

                for (size_t j = 0; j < x_.size(); j++) {
                    indepNew[j] = x_[j];
                }
                Independent(indepNew);

                if (d < newEquations.size() - 1) {
                    // variables with the relationship between x, dxdt and t
                    indep2 = prepareTimeDependentVariables(indepNew);
                } else {
                    indep2 = indepNew;
                    indep2.resize(m);
                }

                Evaluator<Base, CGBase> evaluator(handler0, dep);
                vector<ADCG> depNew = evaluator.evaluate(indep2);

                try {
                    reducedFun_ = new ADFun<CGBase > (indepNew, depNew);
                } catch (const std::exception& ex) {
                    throw CGException(std::string("Failed to create ADFun: ") + ex.what());
                }

#ifdef CPPAD_CG_DAE_VERBOSE
                std::cout << equations.size() << " new equations:\n";
                printModel(reducedFun_);
#endif
            }

            /**
             * Prepare the output information
             */
            newVarInfo.resize(vnodes_.size());
            for (size_t j = 0; j < vnodes_.size(); j++) {
                Vnode<Base>* jj = vnodes_[j];
                size_t tape = jj->tapeIndex();
                if (j < this->varInfo_.size()) {
                    // nothing changed
                    newVarInfo[tape] = this->varInfo_[tape];
                    newVarInfo[tape].setName(jj->name());

                } else {
                    // new variable derivative added by the Pantelides method
                    assert(jj->derivativeOf() != NULL);

                    newVarInfo[tape] = DaeVarInfo(jj->derivativeOf()->tapeIndex(), jj->name());
                }
            }

            std::map<Enode<Base>*, Vnode<Base>*> assigments;
            for (size_t j = 0; j < vnodes_.size(); j++) {
                Vnode<Base>* jj = vnodes_[j];
                if (jj->assigmentEquation() != NULL) {
                    assigments[jj->assigmentEquation()] = jj;
                }
            }

            equationInfo.resize(enodes_.size());
            for (size_t i = 0; i < enodes_.size(); i++) {
                Enode<Base>* ii = enodes_[i];
                int derivativeOf = ii->derivativeOf() != NULL ? ii->derivativeOf()->index() : -1;
                int assignedVarIndex = assigments.count(ii) > 0 ? assigments[ii]->tapeIndex() : -1;

                equationInfo[i] = DaeEquationInfo(i, derivativeOf, assignedVarIndex);
            }
        }

        inline void forwardTimeDiff(const std::vector<Enode<Base>*>& equations,
                                    std::vector<CG<Base> >& dep) const {

            size_t m = reducedFun_->Domain();

            std::vector<CGBase> u(m, CGBase(0));
            u[timeVarIndex_] = CGBase(1);
            std::vector<CGBase> v;
            try {
                v = reducedFun_->Forward(1, u);
            } catch (const std::exception& ex) {
                throw CGException(std::string("Failed to determine model Jacobian (forward mode): ") + ex.what());
            }

            for (size_t e = 0; e < equations.size(); e++) {
                dep[equations[e]->index()] = v[equations[e]->derivativeOf()->index()];
            }
        }

        inline void reverseTimeDiff(const std::vector<Enode<Base>*>& equations,
                                    std::vector<CG<Base> >& dep) const {
            size_t m = reducedFun_->Domain();
            size_t n = reducedFun_->Range();
            std::vector<CGBase> u(m);
            std::vector<CGBase> v(n);

            for (size_t e = 0; e < equations.size(); e++) {
                size_t i = equations[e]->derivativeOf()->index();
                if (reducedFun_->Parameter(i)) { // return zero for this component of f
                    dep[equations[e]->index()] = 0;
                } else {
                    // set v to the i-th coordinate direction
                    v[i] = 1;

                    // compute the derivative of this component of f
                    try {
                        u = reducedFun_->Reverse(1, v);
                    } catch (const std::exception& ex) {
                        throw CGException(std::string("Failed to determine model Jacobian (reverse mode): ") + ex.what());
                    }

                    // reset v to vector of all zeros
                    v[i] = 0;

                    // return the result
                    dep[equations[e]->index()] = u[timeVarIndex_];
                }
            }
        }

        /**
         * Introduces a dependency with respect to time in the provided
         * variables.
         * 
         * \param indepOrig  The variables without time dependency 
         *                    (in the original variable order).
         * \return The new variables with the time dependency 
         *          (in the original variable order).
         */
        inline std::vector<AD<CG<Base> > > prepareTimeDependentVariables(const std::vector<AD<CG<Base> > >& indepOrig) const {
            assert(timeOrigVarIndex_ < vnodes_.size() ||
                   (vnodes_.size() + 1 == indepOrig.size() && timeOrigVarIndex_ == vnodes_.size()));

            using std::vector;
            typedef AD<CGBase> ADCGBase;

            vector<ADCGBase> indepOut(indepOrig.size());
            vector<ADCGBase> ax(3);
            vector<ADCGBase> ay(1);

            assert(timeOrigVarIndex_ < indepOrig.size());
            ax[2] = indepOrig[timeOrigVarIndex_]; // time

            for (size_t j = 0; j < vnodes_.size(); j++) {
                Vnode<Base>* jj = vnodes_[j];
                if (jj->derivative() != NULL) {
                    ax[0] = indepOrig[jj->tapeIndex()]; // x
                    ax[1] = indepOrig[jj->derivative()->tapeIndex()]; // dxdt
                    time_var(0, ax, ay);
                    indepOut[jj->tapeIndex()] = ay[0];
                } else {
                    indepOut[jj->tapeIndex()] = indepOrig[jj->tapeIndex()];
                }
            }

            if (vnodes_.size() < indepOrig.size()) {
                indepOut[indepOut.size() - 1] = indepOrig[timeOrigVarIndex_];
            }

            return indepOut;
        }

        /**
         * Prints out a DAE model to the standard output.
         * 
         * \param fun  The taped model
         */
        inline void printModel(ADFun<CG<Base> >* fun) {
            std::vector<std::string> indepNames(fun->Domain());

            for (size_t j = 0; j < vnodes_.size(); j++) {
                Vnode<Base>* jj = vnodes_[j];
                indepNames[jj->tapeIndex()] = jj->name();
            }

            printModel(fun, indepNames);
        }

        /**
         * Prints out a DAE model to the standard output.
         * 
         * \param fun  The taped model
         * \param vnodes  The independent variables
         */
        inline static void printModel(ADFun<CG<Base> >* fun, const std::vector<std::string>& indepNames) {

            assert(fun != NULL);
            assert(fun->Domain() == indepNames.size() || fun->Domain() == indepNames.size() + 1); // with or without time

            CodeHandler<Base> handler;

            std::vector<CGBase> indep0(fun->Domain());
            handler.makeVariables(indep0);

            std::vector<CGBase> dep0 = fun->Forward(0, indep0);

            CLanguage<double> langC("double");

            /**
             * create variable names
             */
            std::vector<std::string> depNames;

            /**
             * generate the source code
             */
            CLangCustomVariableNameGenerator<double> nameGen(depNames, indepNames,
                                                             "res", "ind", "var");

            std::ostringstream code;
            handler.generateCode(code, langC, dep0, nameGen);
            std::cout << "\n" << code.str() << std::endl;
        }

    private:
        Plantelides(const Plantelides& p); // not implemented
        Plantelides& operator=(const Plantelides& p); // not implemented
    };

}

#endif
