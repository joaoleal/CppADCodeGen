#ifndef CPPAD_CG_PANTELIDES_INCLUDED
#define	CPPAD_CG_PANTELIDES_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad_cgoo/dae_index_reduction/cg_bipartite.hpp>
#include <cppad_cgoo/dae_index_reduction/cg_dae_index_reduction.hpp>
#include <cppad_cgoo/dae_index_reduction/cg_time_diff.hpp>

namespace CppAD {

    /**
     * Pantelides DAE index reduction algorithm
     */
    template<class Base>
    class Plantelides : public DaeIndexReduction<Base> {
    protected:
        // original sparsity
        std::vector<bool> sparsity_;
        // Bipartite graph ([equation i][variable j])
        std::vector<Vnode<Base>*> vnodes_;
        std::vector<Enode<Base>*> enodes_;
        // new index reduced model
        ADFun<CG<Base> >* reducedFun_;
        // the maximum order of the time derivatives in the original model
        size_t origMaxTimeDivOrder_;
    public:

        Plantelides(ADFun<CG<Base> >* fun,
                    const std::vector<DaeVarInfo>& varInfo) :
            DaeIndexReduction<Base>(fun, varInfo),
            reducedFun_(NULL),
            origMaxTimeDivOrder_(0) {

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

            // determine the order of each time derivative

            vector<size_t> derivOrder(n, 0);
            for (size_t dj = 0; dj < n; dj++) {
                int j = varInfo[dj].getDerivativeOf();
                if (j >= 0) {
                    while (j >= 0) {
                        assert(j < n && j != dj);
                        assert(varInfo[j].isTimeDependent());
                        derivOrder[dj]++;
                        j = varInfo[j].getDerivativeOf();
                    }
                }
            }

            origMaxTimeDivOrder_ = *std::max_element(derivOrder.begin(), derivOrder.end());

            // sort the variables according to the time derivative order
            vector<size_t> new2Orig;
            vector<size_t> orig2New(n);
            new2Orig.reserve(n);
            for (size_t order = 0; order <= origMaxTimeDivOrder_; order++) {
                for (size_t j = 0; j < n; j++) {
                    if (derivOrder[j] == order) {
                        orig2New[j] = new2Orig.size();
                        new2Orig.push_back(j);
                    }
                }
            }

            stringstream ss;
            size_t paramCount = 0;
            size_t timeVarCount = 0;
            for (size_t j = 0; j < n; j++) {
                size_t p = new2Orig[j];
                int origIndex = varInfo[p].getDerivativeOf();
                if (origIndex < 0) {
                    // generate the variable name
                    if (varInfo[p].isTimeDependent()) {
                        ss << "x" << timeVarCount;
                        timeVarCount++;
                    } else {
                        ss << "p" << paramCount;
                        paramCount++;
                    }

                    vnodes_[j] = new Vnode<Base > (j, p, ss.str());
                    ss.str("");
                    ss.clear();

                    if (!varInfo[p].isTimeDependent())
                        vnodes_[j]->makeParameter(); // does not depend on time
                } else {
                    vnodes_[j] = new Vnode<Base > (j, p, vnodes_[orig2New[origIndex]]);
                }
            }

            // create the edges
            sparsity_ = jacobianSparsity(*fun);

            for (size_t i = 0; i < m; i++) {
                for (size_t p = 0; p < n; p++) {
                    size_t j = orig2New[p];
                    if (sparsity_[i * n + p] && !vnodes_[j]->isDeleted()) {
                        enodes_[i]->addVariable(vnodes_[j]);
                    }
                }
            }
        }

        /**
         * Performs the DAE differentiation index reductions
         * \param newVarInfo  Variable information of the reduced index model
         * @return the reduced index model
         */
        virtual inline ADFun<CG<Base> >* reduceIndex(std::vector<DaeVarInfo>& newVarInfo) throw (CGException) {
            if (reducedFun_ != NULL)
                throw CGException("reduceIndex() can only be called once!");

            detectSubset2Dif();

#ifdef CPPAD_CG_DAE_VERBOSE
            printResultInfo();
#endif

            generateNewModel(newVarInfo);

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
        inline void generateNewModel(std::vector<DaeVarInfo>& newVarInfo) {
            using std::vector;

            typedef CG<Base> CGBase;
            typedef AD<CGBase> ADCG;

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
                vector<ADCG> indepNew(vnodes_.size() + 1); // variables + time
                Independent(indepNew);

                // variables with the relationship between x dxdt and t
                vector<ADCG> indep2 = prepareTimeDependentVariables(indepNew);
                indep2.resize(indep0.size());

                Evaluator<Base, CG<Base> > evaluator0(handler, dep0);
                vector<ADCG> depNew = evaluator0.evaluate(indep2);
                depNew.resize(enodes_.size());

                reducedFun_ = new ADFun<CGBase > (indepNew, depNew);

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
                size_t n = reducedFun_->Range(); // equation count

                /**
                 * register operations from the other equations
                 */
                CodeHandler<Base> handler0;

                vector<CG<Base> > indep0(m);
                handler0.makeVariables(indep0);

                vector<CG<Base> > dep = reducedFun_->Forward(0, indep0);

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
                    Independent(indepNew);

                    // variables with the relationship between x, dxdt and t
                    indep2 = prepareTimeDependentVariables(indepNew);
                } else {
                    // the very last model creation
                    indepNew.resize(m - 1); // take out time
                    Independent(indepNew);

                    indep2 = indepNew;
                    indep2.resize(m);
                }

                Evaluator<Base, CG<Base> > evaluator(handler0, dep);
                vector<ADCG> depNew = evaluator.evaluate(indep2);

                reducedFun_ = new ADFun<CG<Base> > (indepNew, depNew);

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

                if (jj->derivative() != NULL) {
                    newVarInfo[tape] = jj->derivative()->tapeIndex();
                } else if (jj->derivativeOf() == NULL &&
                        tape < this->varInfo_.size() &&
                        !this->varInfo_[tape].isTimeDependent()) {
                    newVarInfo[tape].makeTimeIndependent();
                } else {
                    newVarInfo[tape] = DaeVarInfo();
                }
            }

        }

        inline void forwardTimeDiff(const std::vector<Enode<Base>*>& equations,
                                    std::vector<CG<Base> >& dep) const {
            typedef CG<Base> CGBase;
            typedef AD<CGBase> ADCG;

            size_t m = reducedFun_->Domain();
            size_t timeIndex = vnodes_.size();

            std::vector<CGBase> u(m, CGBase(0));
            u[timeIndex] = CGBase(1);
            std::vector<CGBase> v = reducedFun_->Forward(1, u);

            for (size_t e = 0; e < equations.size(); e++) {
                dep[equations[e]->index()] = v[equations[e]->derivativeOf()->index()];
            }
        }

        inline void reverseTimeDiff(const std::vector<Enode<Base>*>& equations,
                                    std::vector<CG<Base> >& dep) const {
            size_t m = reducedFun_->Domain();
            size_t n = reducedFun_->Range();
            std::vector<CG<Base> > u(m);
            std::vector<CG<Base> > v(n);

            for (size_t e = 0; e < equations.size(); e++) {
                size_t i = equations[e]->derivativeOf()->index();
                if (reducedFun_->Parameter(i)) { // return zero for this component of f
                    dep[equations[e]->index()] = 0;
                } else {
                    // set v to the i-th coordinate direction
                    v[i] = 1;

                    // compute the derivative of this component of f
                    u = reducedFun_->Reverse(1, v);

                    // reset v to vector of all zeros
                    v[i] = 0;

                    // return the result
                    dep[equations[e]->index()] = u.back();
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
            assert(vnodes_.size() + 1 == indepOrig.size());

            using std::vector;
            typedef AD<CG<Base> > ADCGBase;

            vector<ADCGBase> indepOut(indepOrig.size());
            vector<ADCGBase> ax(3);
            vector<ADCGBase> ay(1);

            ax[2] = indepOrig.back(); // time

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

            indepOut[indepOut.size() - 1] = indepOrig.back(); // time

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

            std::vector<CG<Base> > indep0(fun->Domain());
            handler.makeVariables(indep0);

            std::vector<CG<Base> > dep0 = fun->Forward(0, indep0);

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
