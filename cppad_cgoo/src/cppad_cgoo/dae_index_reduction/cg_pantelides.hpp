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
    public:

        Plantelides(ADFun<CG<Base> >* fun,
                    const std::vector<bool>& eqDifferentialInfo,
                    const std::vector<bool>& varInfo) :
            DaeIndexReduction<Base>(fun, eqDifferentialInfo, varInfo),
            reducedFun_(NULL) {

            // create equation nodes
            enodes_.reserve(1.1 * eqDifferentialInfo.size() + 1);
            enodes_.resize(eqDifferentialInfo.size());
            for (size_t i = 0; i < eqDifferentialInfo.size(); i++) {
                enodes_[i] = new Enode<Base > (i);
            }

            // create the variable nodes
            vnodes_.reserve(2.1 * varInfo.size() + 1);
            vnodes_.resize(varInfo.size());
            std::stringstream ss;
            size_t paramCount = 0;
            size_t timeVarCount = 0;
            for (size_t v = 0; v < varInfo.size(); v++) {
                if (varInfo[v]) {
                    ss << "x" << timeVarCount;
                    timeVarCount++;
                } else {
                    ss << "p" << paramCount;
                    paramCount++;
                }

                vnodes_[v] = new Vnode<Base > (v, ss.str());
                ss.str("");
                ss.clear();

                if (!varInfo[v])
                    vnodes_[v]->makeParameter(); // does not depend on time
            }

            // create the edges
            sparsity_ = jacobianSparsity(*fun);

            size_t m = fun->Range(); // equation count
            size_t n = fun->Domain(); // total variable count

            for (size_t i = 0; i < m; i++) {
                if (eqDifferentialInfo[i]) {
                    assert(varInfo[i]);
                    Vnode<Base>* diffj = new Vnode<Base > (vnodes_.size(), vnodes_[i]);
                    vnodes_.push_back(diffj);
                    enodes_[i]->addVariable(diffj);
                }

                for (size_t j = 0; j < m; j++) {
                    if (sparsity_[i * n + j] && !vnodes_[j]->isDeleted()) {
                        enodes_[i]->addVariable(vnodes_[j]);
                    }
                }
            }
        }

        virtual inline void reduceIndex() {
            detectSubset2Dif();

#ifdef CPPAD_CG_DAE_VERBOSE
            printResultInfo();
#endif

            generateSystem();
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

            delete reducedFun_;
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
                                Vnode<Base>* jDiff = new Vnode<Base > (vnodes_.size(), jj);
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
        inline void dirtyDifferentiateEq(Enode<Base>& i, Enode<Base>& newI) {
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
         * NOT FINISHED!!!!
         */
        inline void generateSystem() {
            using std::vector;

            typedef CG<Base> CGBase;
            typedef AD<CGBase> ADCG;

            vector<vector<Enode<Base>*> > newEquations;

            // find new equations that must be generated by differentiation
            assert(this->eqDifferentialInfo_.size() <= enodes_.size());

            vector<Enode<Base>* > newEqs;
            for (size_t i = 0; i < this->eqDifferentialInfo_.size(); i++) {
                if (enodes_[i]->derivative() != NULL) {
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

            size_t timeIndex = vnodes_.size();

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
                vector<ADCG> indepNew(vnodes_.size() + 1); // variables, time derivatives, time
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

                    // variables with the relationship between x dxdt and t
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
         * Introduces a dependency with respect to time in the provided variables.
         * 
         * \param indepNew  The variables without time dependency
         * \return The new variables with the time dependency
         */
        inline std::vector<AD<CG<Base> > > prepareTimeDependentVariables(const std::vector<AD<CG<Base> > >& indepNew) const {
            assert(vnodes_.size() + 1 == indepNew.size());

            std::vector<AD<CG<Base> > > indep2(indepNew.size());

            std::vector<AD<CG<Base> > > ax(3);
            std::vector<AD<CG<Base> > > ay(1);

            ax[2] = indepNew.back(); // time



            for (size_t j = 0; j < vnodes_.size() - 1; j++) {
                if (vnodes_[j]->derivative() != NULL) {
                    ax[0] = indepNew[j]; // x
                    ax[1] = indepNew[vnodes_[j]->derivative()->index()]; // dxdt
                    time_var(0, ax, ay);
                    indep2[j] = ay[0];
                } else {
                    indep2[j] = indepNew[j];
                }
            }

            indep2[indep2.size() - 1] = indepNew.back(); // time

            return indep2;
        }

        inline std::vector<bool> equationsTimeSparsity(ADFun<CG<Base> >* fun,
                                                       const std::vector<Enode<Base>*>& eqs,
                                                       std::vector<size_t>& row,
                                                       std::vector<size_t>& col) {
            assert(fun != NULL);

            std::vector<bool> sparsity = jacobianForwardSparsity(*fun);

            size_t m = fun->Range(); // equation count
            size_t n = fun->Domain(); // total variable count
            size_t timeIndex = n - 1;

            row.resize(eqs.size());
            col.resize(eqs.size());
            std::fill(col.begin(), col.end(), timeIndex);

#if 1
            std::vector<std::vector<int> > zoot(m, std::vector<int>(n));
            for (size_t i = 0; i < m; i++)
                for (size_t j = 0; j < n; j++)
                    zoot[i][j] = sparsity[i * n + j];
#endif

            for (size_t e = 0; e < eqs.size(); e++) {
                assert(eqs[e]->derivativeOf() != NULL);

                size_t i = eqs[e]->derivativeOf()->index();
                row[e] = i;
                assert(sparsity[i * n + timeIndex]);
            }

            return sparsity;
        }

        /**
         * Prints out a DAE model to the stardard output.
         * 
         * \param fun  The taped model
         */
        inline void printModel(ADFun<CG<Base> >* fun) {

            assert(fun != NULL);

            CodeHandler<Base> handler;

            std::vector<CG<Base> > indep0(fun->Domain());
            handler.makeVariables(indep0);

            std::vector<CG<Base> > dep0 = fun->Forward(0, indep0);

            CLanguage<double> langC("double");

            /**
             * create variable names
             */
            std::vector<std::string> depNames(dep0.size());
            std::vector<std::string> indepNames(indep0.size());

            std::stringstream ss;
            for (size_t i = 0; i < dep0.size(); i++) {
                if (dep0[i].isParameter()) {
                    continue;
                }
                Enode<Base>* ii = enodes_[i];
                size_t depth = 0;
                while (ii->derivativeOf() != NULL) {
                    ii = ii->derivativeOf();
                    depth++;
                }

                if (this->eqDifferentialInfo_[ii->index()]) {
                    for (size_t c = 0; c < depth; c++) {
                        ss << "d";
                    }
                    ss << "dx" << ii->index() << "dt";
                    for (size_t c = 0; c < depth; c++) {
                        ss << "dt";
                    }
                    depNames[i] = ss.str();
                    ss.clear();
                    ss.str("");
                }
            }

            for (size_t j = 0; j < vnodes_.size(); j++) {
                Vnode<Base>* jj = vnodes_[j];
                indepNames[j] = jj->name();
            }

            /**
             * generate the source code
             */
            CLangCustomVariableNameGenerator<double> nameGen(depNames, indepNames,
                                                             "res", "ind", "var");

            std::ostringstream code;
            handler.generateCode(code, langC, dep0, nameGen);
            std::cout << "\n" << code.str() << std::endl;
        }

    };

}

#endif

