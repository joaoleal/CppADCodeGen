#ifndef CPPAD_CG_LOOP_INCLUDED
#define CPPAD_CG_LOOP_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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

namespace CppAD {

    template<class Base>
    class IndependentOrder {
    public:
        // provides an independent variable for each loop iteration
        const std::vector<const OperationNode<Base>*> order;

        inline IndependentOrder(const std::vector<const OperationNode<Base>*>& myOrder) :
            order(myOrder) {
        }
    };

    template<class Base>
    class OperationArgumentsIndepOrder {
    public:
        std::map<size_t, IndependentOrder<Base>*> arg2Order;
    };

    /**
     * 
     */
    template<class Base>
    class LoopCodeHandler : public CodeHandler<Base> {
    public:
        using CodeHandler<Base>::manageOperationNode;

    };

    /**
     * A for loop.
     */
    template<class Base>
    class Loop {
    private:
        typedef std::pair<size_t, CG<Base> > IndexValue;
    public:
        /**
         * The equations inside the loop
         */
        std::set<EquationPattern<Base>*> equations;
        /**
         * Which argument positions of operations (from the reference dependent)
         * use indexed independent variables 
         * (operation -> argument index -> iteration -> independent)
         */
        IndexedIndependent<Base> indexedOpIndep;
    private:
        /**
         * The number of iterations this loop will have
         */
        size_t iterationCount_;
        /**
         * code handler for the operations in the loop
         */
        LoopCodeHandler<Base> handler_;
        /**
         * Groups of equations which share common temporary variables
         */
        std::vector<EquationGroup<Base> > eqGroups_;
        /**
         * The order of equation patterns in the loop's tape
         */
        std::map<EquationPattern<Base>*, size_t> equationOrder_;
        /**
         * The evaluated dependents in each loop iteration
         * ([iteration 1]{dep1, dep3, ...}; [iteration 2]{dep5, dep6, ...}; ...)
         */
        std::vector<std::set<size_t> > iterationDependents_;
        /**
         * Map each dependent to its corresponding iteration
         */
        std::map<size_t, size_t> dep2Iteration_;
        /**
         * total number of independent variables in the loop's tape
         */
        size_t nIndependents_;
        /**
         * 
         */
        LoopModel<Base>* loopModel_;
        /**
         * indexed independent variables for the loop (clones)
         */
        std::map<const OperationNode<Base>*, IndexValue> independentsIndexed_;
        /**
         * non-indexed independent variables for the loop (clones)
         */
        std::map<const OperationNode<Base>*, IndexValue> independentsNonIndexed_;
        /**
         * independent variables for the loop created from temporary variables (clones)
         */
        std::map<const OperationNode<Base>*, IndexValue> independentsTemp_;
        /**
         * id -> clone
         */
        std::map<size_t, OperationNode<Base>*> clonesTemporary_;
        /**
         * orig -> clone
         */
        std::map<const OperationNode<Base>*, OperationNode<Base>*> temporaryClone2Orig_;
        /**
         * orig -> clone
         */
        std::map<const OperationNode<Base>*, OperationNode<Base>*> orig2ConstIndepClone_;
        /**
         * independent order -> clone
         */
        std::map<const IndependentOrder<Base>*, OperationNode<Base>*> indexedIndep2clone_;
        /**
         * Maps the independent from the first loop iteration to its independent
         * variable order
         * (first independent in the order -> all possible orders)
         */
        std::map<const OperationNode<Base>*, std::vector<IndependentOrder<Base>*> > firstIndep2orders_;
        /**
         * Maps an operations to their arguments and the corresponding
         * independent variable order
         */
        std::map<const OperationNode<Base>*, OperationArgumentsIndepOrder<Base>* > op2Arg2IndepOrder_;
        /**
         * used just to keep pointers so that the objects can be deleted later
         */
        std::list<OperationArgumentsIndepOrder<Base>*> arg2IndepOrder_; // consider changing to forward_list (C++11)
        /**
         * Variable id counter
         */
        size_t idCounter_;
    public:

        inline Loop(EquationPattern<Base>& eq) :
            iterationCount_(0), // not known yet (only after all equations have been added)
            eqGroups_(1),
            nIndependents_(0),
            loopModel_(NULL) {
            //indexedOpIndep(eq.indexedOpIndep), // must be determined later for a different reference
            //constOperationIndependents(eq.constOperationIndependents) {
            equations.insert(&eq);
            eqGroups_[0].equations.insert(&eq);
        }

        inline void addEquation(EquationPattern<Base>& eq) {
            equations.insert(&eq);
            eqGroups_[0].equations.insert(&eq);
        }

        inline void setLinkedDependents(const std::set<std::set<size_t>*>& newLoopRelations) {
            assert(eqGroups_.size() == 1);

            eqGroups_[0].linkedDependents.clear();
            eqGroups_[0].linkedDependents.reserve(newLoopRelations.size());

            std::set<std::set<size_t>*>::const_iterator it;
            for (it = newLoopRelations.begin(); it != newLoopRelations.end(); ++it)
                eqGroups_[0].linkedDependents.push_back(**it);
        }

        inline void addLinkedEquationsByNonIndexed(EquationPattern<Base>* eq1,
                                                   EquationPattern<Base>* eq2) {
            eqGroups_[0].addLinkedEquationsByNonIndexed(eq1, eq2);
        }

        inline size_t getLinkedEquationsByNonIndexedCount() const {
            return eqGroups_[0].getLinkedEquationsByNonIndexedCount();
        }

        /**
         * The number of iterations this loop will have
         */
        inline size_t getIterationCount() const {
            return iterationCount_;
        }

        const std::vector<std::set<size_t> >& getIterationDependents() const {
            return iterationDependents_;
        }

        inline LoopModel<Base>* getModel() const {
            return loopModel_;
        }

        inline LoopModel<Base>* releaseLoopModel() {
            LoopModel<Base>* loopAtomic = loopModel_;
            loopModel_ = NULL;
            return loopAtomic;
        }

        /**
         * Combines the provided loop with the current one
         * 
         * @param other The other loop
         */
        void merge(Loop<Base>& other,
                   const std::set<EquationPattern<Base>*>& indexedLoopRelations,
                   const std::vector<std::pair<EquationPattern<Base>*, EquationPattern<Base>*> >& nonIndexedLoopRelations) {
            assert(iterationCount_ == other.iterationCount_);


            equations.insert(other.equations.begin(), other.equations.end());
            other.equations.clear(); // so that it does not delete the equations

            assert(eqGroups_.size() == 1);
            assert(other.eqGroups_.size() == 1);

            EquationGroup<Base>& g = eqGroups_[0];
            EquationGroup<Base>& og = other.eqGroups_[0];
            g.equations.insert(og.equations.begin(), og.equations.end());

            assert(equationOrder_.empty());
            assert(iterationDependents_.empty());

            g.linkedEquationsByNonIndexed.insert(og.linkedEquationsByNonIndexed.begin(), og.linkedEquationsByNonIndexed.end());
            typename std::set<EquationPattern<Base>*>::const_iterator itNIndexed;
            for (itNIndexed = indexedLoopRelations.begin(); itNIndexed != indexedLoopRelations.end(); ++itNIndexed) {
                g.linkedEquationsByNonIndexed.erase(*itNIndexed);
            }

            for (size_t e = 0; e < nonIndexedLoopRelations.size(); e++) {
                g.addLinkedEquationsByNonIndexed(nonIndexedLoopRelations[e].first, nonIndexedLoopRelations[e].second);
            }

        }

        void mergeEqGroups(Loop<Base>& other) {
            eqGroups_.insert(eqGroups_.end(), other.eqGroups_.begin(), other.eqGroups_.end());

            size_t nEq = equations.size();
            equations.insert(other.equations.begin(), other.equations.end());
            other.equations.clear(); // so that it does not delete the equations

            /**
             * Update equation index
             */
            typename std::map<EquationPattern<Base>*, size_t>::const_iterator it;
            for (it = other.equationOrder_.begin(); it != other.equationOrder_.end(); ++it) {
                equationOrder_[it->first] = it->second + nEq;
            }

            assert(iterationDependents_.size() == other.iterationDependents_.size());
            for (size_t iter = 0; iter < iterationDependents_.size(); iter++) {
                iterationDependents_[iter].insert(other.iterationDependents_[iter].begin(), other.iterationDependents_[iter].end());
            }

        }

        void createLoopModel(const std::vector<CG<Base> >& dependents,
                             const std::vector<CG<Base> >& independents,
                             const std::map<size_t, EquationPattern<Base>*>& dep2Equation,
                             std::map<OperationNode<Base>*, size_t>& origTemp2Index) throw (CGException) {

            assert(dep2Iteration_.empty());
            for (size_t iter = 0; iter < iterationCount_; iter++) {
                const std::set<size_t>& deps = iterationDependents_[iter];

                std::set<size_t>::const_iterator itDeps;
                for (itDeps = deps.begin(); itDeps != deps.end(); ++itDeps) {
                    dep2Iteration_[*itDeps] = iter;
                }
            }

            /**
             * Determine the reference iteration for each
             */
            for (size_t g = 0; g < eqGroups_.size(); g++) {
                eqGroups_[g].findReferenceIteration();
#ifdef CPPADCG_PRINT_DEBUG
                std::cout << "reference iteration=" << eqGroups_[g].refIteration << "\n";
                print(eqGroups_[g].iterationDependents);
                std::cout << std::endl;

                typename std::set<EquationPattern<Base>*>::const_iterator iteq;
                for (iteq = eqGroups_[g].equations.begin(); iteq != eqGroups_[g].equations.end(); ++iteq) {
                    std::cout << "eq dependents=";
                    print((*iteq)->dependents);
                    std::cout << std::endl;
                }
#endif
            }


            std::set<size_t>::const_iterator itDep;

            for (size_t g = 0; g < eqGroups_.size(); g++) {
                const EquationGroup<Base>& group = eqGroups_[g];
                const std::set<size_t>& refItDep = group.iterationDependents[group.refIteration];
                assert(refItDep.size() == group.equations.size());

                for (itDep = refItDep.begin(); itDep != refItDep.end(); ++itDep) {
                    size_t dep = *itDep;
                    EquationPattern<Base>::uncolor(dependents[dep].getOperationNode());
                }

                for (itDep = refItDep.begin(); itDep != refItDep.end(); ++itDep) {
                    size_t dep = *itDep;
                    EquationPattern<Base>* eq = dep2Equation.at(dep);

                    // operations that use indexed independent variables in the reference iteration
                    std::set<const OperationNode<Base>*> indexedOperations;

                    eq->colorIndexedPath(dep, dependents, 1, indexedOperations);
                    if (dep == eq->depRefIndex) {
                        const std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >& op2Arguments = eq->indexedOpIndep.op2Arguments;
                        typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::const_iterator itop2a;
                        for (itop2a = op2Arguments.begin(); itop2a != op2Arguments.end(); ++itop2a) {
                            // currently there is no way to make a distinction between yi = xi and y_(i+1) = x_(i+1)
                            // since both operations which use indexed independents would be NULL (the dependent)
                            // an alias is used for these cases
                            assert(itop2a->first != NULL);
                            addOperationArguments2Loop(itop2a->first, itop2a->second);
                        }

                    } else {
                        // generate loop references
                        typename std::set<const OperationNode<Base>*>::const_iterator it;
                        for (it = indexedOperations.begin(); it != indexedOperations.end(); ++it) {
                            const OperationNode<Base>* opLoopRef = *it;
                            // currently there is no way to make a distinction between yi = xi and y_(i+1) = x_(i+1)
                            // since both operations which use indexed independents would be NULL (the dependent)
                            // an alias is used for these cases
                            assert(opLoopRef != NULL);

                            const OperationNode<Base>* opEqRef = eq->operationEO2Reference.at(dep).at(opLoopRef); /////////////////////////////////////////////////////
                            addOperationArguments2Loop(opLoopRef, eq->indexedOpIndep.op2Arguments.at(opEqRef));
                        }
                    }

                    // not need anymore (lets free this memory now)
                    eq->indexedOpIndep.op2Arguments.clear();
                }
            }

            /**
             * independent variable index patterns
             */
            generateIndependentLoopIndexes();

            resetCounters(dependents);

            createLoopTapeNModel(dependents, independents, dep2Equation, origTemp2Index);

            /**
             * Clean-up
             */
            for (size_t g = 0; g < eqGroups_.size(); g++) {
                const EquationGroup<Base>& group = eqGroups_[g];
                const std::set<size_t>& refItDep = group.iterationDependents[group.refIteration];
                for (itDep = refItDep.begin(); itDep != refItDep.end(); ++itDep) {
                    size_t dep = *itDep;
                    EquationPattern<Base>::uncolor(dependents[dep].getOperationNode());
                }
            }
            resetCounters(dependents);
        }

        void generateDependentLoopIndexes(const std::map<size_t, EquationPattern<Base>*>& dep2Equation) {
            iterationDependents_.clear();
            equationOrder_.clear();
            iterationCount_ = 0;
            size_t nMaxIt = 0;
            /**
             * assign a dependent variable from each equation to an iteration
             */
            std::map<EquationPattern<Base>*, std::set<size_t> > depsInEq;

            typename std::set<EquationPattern<Base>*>::const_iterator eqIt;
            for (eqIt = equations.begin(); eqIt != equations.end(); ++eqIt) {
                EquationPattern<Base>* eq = *eqIt;
                depsInEq[eq] = eq->dependents;
                nMaxIt = std::max(nMaxIt, eq->dependents.size());
            }

            iterationDependents_.reserve(nMaxIt + 2 * equations.size());

            for (size_t g = 0; g < eqGroups_.size(); g++) {
                EquationGroup<Base>& group = eqGroups_[g];
                const std::set<EquationPattern<Base>*>& eqs = group.equations;
                const std::vector<std::set<size_t> >& linkedDependents = group.linkedDependents;
                std::vector<std::set<size_t> >& relatedEqIterationDeps = group.iterationDependents;

                relatedEqIterationDeps.reserve(iterationDependents_.size());

                // assign an index to each equation
                typename std::set<EquationPattern<Base>*>::const_iterator eqIt;
                for (eqIt = eqs.begin(); eqIt != eqs.end(); ++eqIt) {
                    EquationPattern<Base>* eq = *eqIt;
                    size_t eqo_size = equationOrder_.size();
                    equationOrder_[eq] = eqo_size;
                }

                // sort dependents
                std::set<size_t> dependents;
                for (eqIt = eqs.begin(); eqIt != eqs.end(); ++eqIt) {
                    EquationPattern<Base>* eq = *eqIt;
                    dependents.insert(eq->dependents.begin(), eq->dependents.end());
                }

                std::map<size_t, std::map<EquationPattern<Base>*, std::set<size_t> > > nIndexedGroupPos2Eq2deps;

                std::set<size_t>::const_iterator itDep = dependents.begin();
                size_t i = 0; // iteration counter
                while (itDep != dependents.end()) {
                    size_t dep = *itDep;
                    itDep++;

                    if (iterationDependents_.size() <= i)
                        iterationDependents_.resize(i + 1);
                    std::set<size_t>& itDepi = iterationDependents_[i];

                    if (relatedEqIterationDeps.size() <= i)
                        relatedEqIterationDeps.resize(i + 1);
                    std::set<size_t>& ritDepi = relatedEqIterationDeps[i];

                    long pos = group.findIndexedLinkedDependent(dep);
                    if (pos >= 0) {
                        // assign the dependent to the first iteration with all its relationships
                        std::set<size_t>::const_iterator itDep2;
                        for (itDep2 = linkedDependents[pos].begin(); itDep2 != linkedDependents[pos].end(); ++itDep2) {
                            size_t dep2 = *itDep2;
                            if (dep2 == *itDep) itDep++; //make sure the iterator is valid
                            itDepi.insert(dep2);
                            ritDepi.insert(dep2);
                            dependents.erase(dep2);
                        }

                        i++;
                    } else {
                        // maybe this dependent shares a non-indexed temporary variable
                        EquationPattern<Base>* eq = dep2Equation.at(dep);

                        long posN = group.findNonIndexedLinkedRel(eq);
                        if (posN >= 0) {
                            // there are only non-indexed shared relations with other equations (delay processing...)
                            dependents.erase(dep);
                            nIndexedGroupPos2Eq2deps[posN][eq].insert(dep);
                        } else {
                            itDepi.insert(dep);
                            ritDepi.insert(dep);

                            i++;
                        }
                    }

                }

                /**
                 * place dependents which only share non-indexed variables
                 */
                if (!nIndexedGroupPos2Eq2deps.empty()) {

                    std::map<EquationPattern<Base>*, std::set<size_t> > eqIterations;
                    for (size_t i = 0; i < relatedEqIterationDeps.size(); i++) {
                        const std::set<size_t>& deps = relatedEqIterationDeps[i];
                        std::set<size_t>::const_iterator itDep;
                        for (itDep = deps.begin(); itDep != deps.end(); ++itDep) {
                            size_t dep = *itDep;
                            eqIterations[dep2Equation.at(dep)].insert(i);
                        }
                    }

                    typename std::map<size_t, std::map<EquationPattern<Base>*, std::set<size_t> > >::iterator itPos2Eq2Dep;
                    for (itPos2Eq2Dep = nIndexedGroupPos2Eq2deps.begin(); itPos2Eq2Dep != nIndexedGroupPos2Eq2deps.end(); ++itPos2Eq2Dep) {
                        size_t posN = itPos2Eq2Dep->first;
                        // must pick one dependent from each equation for each iteration
                        std::vector<size_t> deps;
                        deps.reserve(itPos2Eq2Dep->second.size());

                        std::set<size_t> usedIterations; // iterations used by these equations 
                        // determine used iteration indexes
                        const std::set<EquationPattern<Base>*>& relations = group.linkedEquationsByNonIndexedRel[posN];
                        typename std::set<EquationPattern<Base>*>::const_iterator itRel;
                        for (itRel = relations.begin(); itRel != relations.end(); ++itRel) {
                            const std::set<size_t>& iters = eqIterations[*itRel];
                            usedIterations.insert(iters.begin(), iters.end());
                        }
                        
                        while (true) {
                            
                            // must pick one dependent from each equation for each iteration
                            deps.clear();
                            
                            typename std::map<EquationPattern<Base>*, std::set<size_t> >::iterator itEq2Dep;
                            for (itEq2Dep = itPos2Eq2Dep->second.begin(); itEq2Dep != itPos2Eq2Dep->second.end(); ++itEq2Dep) {
                                if (!itEq2Dep->second.empty()) {
                                    deps.push_back(*itEq2Dep->second.begin());
                                    itEq2Dep->second.erase(itEq2Dep->second.begin());
                                }
                            }

                            if (deps.empty()) {
                                break; // done
                            }

                            // find a free iteration index
                            size_t i = 0;
                            std::set<size_t>::const_iterator itIter;
                            for (itIter = usedIterations.begin(); itIter != usedIterations.end();) {
                                size_t i1 = *itIter;
                                ++itIter;
                                if (itIter != usedIterations.end()) {
                                    size_t i2 = *itIter;
                                    if (i2 - i1 != 1) {
                                        i = i1 + 1;
                                        break;
                                    }
                                } else {
                                    i = i1 + 1;
                                }
                            }

                            // add the dependents to the iteration
                            usedIterations.insert(i);

                            if (iterationDependents_.size() <= i)
                                iterationDependents_.resize(i + 1);
                            std::set<size_t>& itDepi = iterationDependents_[i];

                            if (relatedEqIterationDeps.size() <= i)
                                relatedEqIterationDeps.resize(i + 1);
                            std::set<size_t>& ritDepi = relatedEqIterationDeps[i];
                            itDepi.insert(deps.begin(), deps.end());
                            ritDepi.insert(deps.begin(), deps.end());
                        }
                    }
                    /**
                     * @todo reorder iterations according to the lowest
                     *       dependent in each iteration if there were new
                     *       iterations only with dependents related by
                     *       non-indexed shared variables
                     */

                }
                
                iterationCount_ = std::max(iterationCount_, iterationDependents_.size());
            }

        }

        /**
         * Destructor
         */
        virtual ~Loop() {
            typename std::map<const OperationNode<Base>*, std::vector<IndependentOrder<Base>*> >::const_iterator itf;
            for (itf = firstIndep2orders_.begin(); itf != firstIndep2orders_.end(); ++itf) {
                const std::vector<IndependentOrder<Base>*>& orders = itf->second;
                typename std::vector<IndependentOrder<Base>*>::const_iterator ito;
                for (ito = orders.begin(); ito != orders.end(); ++ito) {
                    delete *ito;
                }
            }

            typename std::list<OperationArgumentsIndepOrder<Base>*>::const_iterator itio;
            for (itio = arg2IndepOrder_.begin(); itio != arg2IndepOrder_.end(); ++itio) {
                delete *itio;
            }

            typename std::set<EquationPattern<Base>*>::const_iterator it;
            for (it = equations.begin(); it != equations.end(); ++it) {
                delete *it;
            }

            delete loopModel_;
        }

        /***********************************************************************
         *                            private
         **********************************************************************/
    private:

        void addOperationArguments2Loop(const OperationNode<Base>* op,
                                        const OperationIndexedIndependents<Base>& eqOpIndeIndep) {
            assert(!dep2Iteration_.empty());

            OperationIndexedIndependents<Base>& loopOpIndeIndep = indexedOpIndep.op2Arguments[op];
            loopOpIndeIndep.arg2Independents.resize(eqOpIndeIndep.arg2Independents.size());

            // some iterations might have not been defined by previous equation patterns
            for (size_t a = 0; a < eqOpIndeIndep.arg2Independents.size(); a++) {
                if (eqOpIndeIndep.arg2Independents[a].empty())
                    continue;

                typename std::map<size_t, const OperationNode<Base>*>::const_iterator itDepIndep;
                for (itDepIndep = eqOpIndeIndep.arg2Independents[a].begin(); itDepIndep != eqOpIndeIndep.arg2Independents[a].end(); ++itDepIndep) {
                    size_t dep = itDepIndep->first;
                    size_t iter = dep2Iteration_.at(dep);
                    loopOpIndeIndep.arg2Independents[a][iter] = itDepIndep->second;
                }
            }

        }

        void generateIndependentLoopIndexes() {
            assert(iterationCount_ > 0); //number of iterations and dependent indexes must have already been determined

            // loop all operations from the reference dependents which use indexed independents
            typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::const_iterator it;
            for (it = indexedOpIndep.op2Arguments.begin(); it != indexedOpIndep.op2Arguments.end(); ++it) {
                const OperationNode<Base>* operation = it->first;
                const OperationIndexedIndependents<Base>& opInd = it->second;

                OperationArgumentsIndepOrder<Base>* arg2orderPos = new OperationArgumentsIndepOrder<Base>();
                op2Arg2IndepOrder_[operation] = arg2orderPos;

                arg2IndepOrder_.push_front(arg2orderPos);

                // loop all arguments
                size_t aSize = opInd.arg2Independents.size();
                for (size_t argumentIndex = 0; argumentIndex < aSize; argumentIndex++) {
                    const std::map<size_t, const OperationNode<Base>*>& dep2Indep = opInd.arg2Independents[argumentIndex];
                    if (dep2Indep.empty())
                        continue; // not an indexed variable

                    std::vector<const OperationNode<Base>*> order(iterationCount_);

                    /**
                     * create the independent variable order
                     */
                    assert(dep2Indep.size() > 0 && dep2Indep.size() <= iterationCount_);
                    typename std::map<size_t, const OperationNode<Base>*>::const_iterator itDep2Indep;
                    for (itDep2Indep = dep2Indep.begin(); itDep2Indep != dep2Indep.end(); ++itDep2Indep) {
                        size_t iterationIndex = itDep2Indep->first;
                        const OperationNode<Base>* indep = itDep2Indep->second;

                        order[iterationIndex] = indep;
                    }

                    /**
                     * try to find an existing independent variable order
                     */
                    std::vector<IndependentOrder<Base>*>& availableOrders = firstIndep2orders_[order[0]];
                    IndependentOrder<Base>* match = NULL;

                    long a_size = availableOrders.size();
                    for (long o = 0; o < a_size; o++) {
                        IndependentOrder<Base>* orderO = availableOrders[o];
                        bool ok = true;
                        for (size_t iterationIndex = 0; iterationIndex < iterationCount_; iterationIndex++) {
                            if (orderO->order[iterationIndex] != order[iterationIndex]) {
                                ok = false;
                                break;
                            }
                        }
                        if (ok) {
                            match = orderO;
                            break;
                        }
                    }

                    if (match != NULL) {
                        // found another operation with the same independent variable order
                        arg2orderPos->arg2Order[argumentIndex] = match;
                    } else {
                        // brand new independent variable order
                        IndependentOrder<Base>* iOrder = new IndependentOrder<Base>(order);
                        availableOrders.push_back(iOrder);
                        arg2orderPos->arg2Order[argumentIndex] = iOrder;
                    }

                }

            }
        }

        /**
         * Resets the node counters for the reference iteration of each equation
         * group
         * 
         * @param dependents
         */
        void resetCounters(const std::vector<CG<Base> >& dependents) {
            for (size_t g = 0; g < eqGroups_.size(); g++) {
                const EquationGroup<Base>& group = eqGroups_[g];
                const std::set<size_t>& itDeps = group.iterationDependents[group.refIteration];

                std::set<size_t>::const_iterator depIt;
                for (depIt = itDeps.begin(); depIt != itDeps.end(); ++depIt) {
                    size_t depIndex = *depIt;
                    OperationNode<Base>* node = dependents[depIndex].getOperationNode();
                    if (node != NULL) {
                        Loop<Base>::resetCounters(*node);
                    }
                }
            }
        }

        /**
         * Creates the model for the loop equations
         * 
         * @param dependents original model dependent variable vector
         * @param independents original model independent variable vector
         * @param dep2Equation maps an equation/dependent index to an equation pattern
         * @param origTemp2Index 
         */
        void createLoopTapeNModel(const std::vector<CG<Base> >& dependents,
                                  const std::vector<CG<Base> >& independents,
                                  const std::map<size_t, EquationPattern<Base>*>& dep2Equation,
                                  std::map<OperationNode<Base>*, size_t>& origTemp2Index) throw (CGException) {
            typedef CG<Base> CGB;
            typedef AD<CGB> ADCGB;
            assert(independents.size() > 0);
            assert(independents[0].getCodeHandler() != NULL);
            CodeHandler<Base>& origHandler = *independents[0].getCodeHandler();

            /**
             * create the new/clone operations for the reference iteration only
             */
            nIndependents_ = 0;
            idCounter_ = 1;

            assert(equationOrder_.size() == equations.size());

            std::vector<CGB> deps(equations.size());

            for (size_t g = 0; g < eqGroups_.size(); g++) {
                const EquationGroup<Base>& group = eqGroups_[g];
                const std::set<size_t>& iterationDependents = group.iterationDependents[group.refIteration];

                std::set<size_t>::const_iterator depIt;
                for (depIt = iterationDependents.begin(); depIt != iterationDependents.end(); ++depIt) {
                    size_t depIndex = *depIt;
                    EquationPattern<Base>* eq = dep2Equation.at(depIndex);
                    OperationNode<Base>* node = dependents[depIndex].getOperationNode();

                    Argument<Base> aClone;

                    if (node != NULL) {
                        if (node->getOperationType() == CGInvOp) {
                            aClone = Argument<Base>(createIndependentClone(NULL, 0, *node));
                        } else {
                            aClone = makeGraphClones(*eq, *node);
                        }
                    } else {
                        aClone = Argument<Base>(dependents[depIndex].getValue());
                    }

                    size_t i = equationOrder_.at(eq);
                    deps[i] = CGB(handler_, aClone);
                }
            }

            /*******************************************************************
             * create the tape for the reference iteration
             ******************************************************************/
            std::map<const OperationNode<Base>*, size_t> origModelIndepOrder;
            size_t nOrigIndep = independents.size();
            for (size_t j = 0; j < nOrigIndep; j++) {
                origModelIndepOrder[independents[j].getOperationNode()] = j;
            }

            // indexed independents
            assert(indexedIndep2clone_.size() == independentsIndexed_.size());

            std::vector<const OperationNode<Base>*> indexedCloneOrder;
            indexedCloneOrder.reserve(indexedIndep2clone_.size());

            std::map<const OperationNode<Base>*, const IndependentOrder<Base>*> clone2indexedIndep;
            typename std::map<const IndependentOrder<Base>*, OperationNode<Base>*>::const_iterator it;
            for (it = indexedIndep2clone_.begin(); it != indexedIndep2clone_.end(); ++it) {
                clone2indexedIndep[it->second] = it->first;
                indexedCloneOrder.push_back(it->second);
            }

            struct IndexedIndepSorter indexedSorter(clone2indexedIndep, origModelIndepOrder);
            std::sort(indexedCloneOrder.begin(), indexedCloneOrder.end(), indexedSorter);

            // original indep index -> non-indexed independent clones
            std::map<size_t, const OperationNode<Base>*> nonIndexedCloneOrder;

            // [tape variable] = original independent index
            std::map<const OperationNode<Base>*, const OperationNode<Base>*> clones2ConstIndep;
            typename std::map<const OperationNode<Base>*, OperationNode<Base>*>::const_iterator itc;
            size_t s = 0;
            for (itc = orig2ConstIndepClone_.begin(); itc != orig2ConstIndepClone_.end(); ++itc, s++) {
                const OperationNode<Base>* orig = itc->first;
                OperationNode<Base>* clone = itc->second;

                size_t j = origModelIndepOrder.at(orig);
                clones2ConstIndep[clone] = orig;
                nonIndexedCloneOrder[j] = clone;
            }

            size_t nIndexed = indexedCloneOrder.size();
            size_t nNonIndexed = nonIndexedCloneOrder.size();
            size_t nTmpIndexed = independentsTemp_.size();

            // tape independent array
            size_t nIndep = independentsIndexed_.size() +
                    independentsNonIndexed_.size() +
                    independentsTemp_.size();
            std::vector<ADCGB> loopIndeps(nIndep);

            typename std::map<const OperationNode<Base>*, IndexValue>::const_iterator itt;
            typename std::map<size_t, const OperationNode<Base>*>::const_iterator origJ2CloneIt;

            if (nIndep == 0) {
                loopIndeps.resize(1); // the tape cannot have 0 independents
                loopIndeps[0] = Base(0);

            } else {

                /**
                 * assign values from the original model if possible (to avoid NaN)
                 */
                for (size_t j = 0; j < nIndexed; j++) {
                    const IndexValue& iv = independentsIndexed_.at(indexedCloneOrder[j]);
                    loopIndeps[j] = iv.second;
                }

                s = 0;
                for (origJ2CloneIt = nonIndexedCloneOrder.begin(); origJ2CloneIt != nonIndexedCloneOrder.end(); ++origJ2CloneIt, s++) {
                    const IndexValue& iv = independentsNonIndexed_.at(origJ2CloneIt->second);
                    loopIndeps[nIndexed + s] = iv.second;
                }

                s = 0;
                for (itt = independentsTemp_.begin(); itt != independentsTemp_.end(); ++itt, s++) {
                    const IndexValue& iv = itt->second;
                    loopIndeps[nIndexed + nNonIndexed + s] = iv.second;
                }
            }

            /**
             * register the independent variables
             */
            CppAD::Independent(loopIndeps);

            /**
             * Reorder independent variables for the new tape (reference iteration only)
             */
            std::vector<ADCGB> localIndeps(nIndep);
            if (nIndep > 0) {
                for (size_t j = 0; j < nIndexed; j++) {
                    const IndexValue& iv = independentsIndexed_.at(indexedCloneOrder[j]);
                    size_t localIndex = iv.first;
                    localIndeps[localIndex] = loopIndeps[j];
                }

                s = 0;
                for (origJ2CloneIt = nonIndexedCloneOrder.begin(); origJ2CloneIt != nonIndexedCloneOrder.end(); ++origJ2CloneIt, s++) {
                    size_t localIndex = independentsNonIndexed_.at(origJ2CloneIt->second).first;
                    localIndeps[localIndex] = loopIndeps[nIndexed + s];
                }

                s = 0;
                for (itt = independentsTemp_.begin(); itt != independentsTemp_.end(); ++itt, s++) {
                    size_t localIndex = itt->second.first;
                    localIndeps[localIndex] = loopIndeps[nIndexed + nNonIndexed + s];
                }
            }

            /**
             * create the tape
             */
            Evaluator<Base, CGB> evaluator1stIt(handler_, deps);

            // load any atomic function used in the original model
            const std::map<size_t, CGAbstractAtomicFun<Base>* >& atomicsOrig = origHandler.getAtomicFunctions();
            std::map<size_t, atomic_base<CGB>* > atomics;
            atomics.insert(atomicsOrig.begin(), atomicsOrig.end());
            evaluator1stIt.addAtomicFunctions(atomics);

            std::vector<ADCGB> newDeps = evaluator1stIt.evaluate(localIndeps);

            std::auto_ptr<ADFun<CGB> >funIndexed(new ADFun<CGB>());
            funIndexed->Dependent(newDeps);

            /*******************************************************************
             * create the atomic loop function
             ******************************************************************/
            std::vector<std::vector<size_t> > dependentOrigIndexes(equations.size(), std::vector<size_t> (iterationCount_, std::numeric_limits<size_t>::max()));
            for (size_t it = 0; it < iterationCount_; it++) {
                const std::set<size_t>& itDeps = iterationDependents_[it];
                std::set<size_t>::const_iterator itDep;
                for (itDep = itDeps.begin(); itDep != itDeps.end(); ++itDep) {
                    size_t origDep = *itDep;
                    EquationPattern<Base>* eq = dep2Equation.at(origDep);
                    size_t i = equationOrder_.at(eq);
                    dependentOrigIndexes[i][it] = origDep;
                }
            }

            //[tape variable][iteration] =  original independent index
            std::vector<std::vector<size_t> > indexedIndependents(nIndexed, std::vector<size_t>(iterationCount_));
            for (size_t j = 0; j < indexedCloneOrder.size(); j++) {
                const IndependentOrder<Base>* origOrder = clone2indexedIndep.at(indexedCloneOrder[j]);
                for (size_t it = 0; it < iterationCount_; it++) {
                    const OperationNode<Base>* indep = origOrder->order[it];
                    size_t index;
                    if (indep != NULL) {
                        index = origModelIndepOrder.at(indep);
                    } else {
                        index = std::numeric_limits<size_t>::max(); // not used at this iteration by any equation
                    }
                    indexedIndependents[j][it] = index;
                }
            }

            std::vector<size_t> temporaryIndependents(nTmpIndexed);

            size_t j = 0;
            for (itt = independentsTemp_.begin(); itt != independentsTemp_.end(); ++itt, j++) {
                const OperationNode<Base>* tmpClone = itt->first;
                OperationNode<Base>* origTmpNode = temporaryClone2Orig_.at(tmpClone);

                /**
                 * assign an index (k) to each temporary variable 
                 */
                size_t k;
                typename std::map<OperationNode<Base>*, size_t>::const_iterator itz = origTemp2Index.find(origTmpNode);
                if (itz == origTemp2Index.end()) {
                    k = origTemp2Index.size();
                    origTemp2Index[origTmpNode] = k; // new index for a temporary variable
                } else {
                    k = itz->second;
                }

                temporaryIndependents[j] = k;
            }

            std::vector<size_t> nonIndexedIndependents(orig2ConstIndepClone_.size());
            s = 0;
            for (origJ2CloneIt = nonIndexedCloneOrder.begin(); origJ2CloneIt != nonIndexedCloneOrder.end(); ++origJ2CloneIt, s++) {
                nonIndexedIndependents[s] = origJ2CloneIt->first;
            }

            loopModel_ = new LoopModel<Base>(funIndexed.release(),
                    iterationCount_,
                    dependentOrigIndexes,
                    indexedIndependents,
                    nonIndexedIndependents,
                    temporaryIndependents);

            loopModel_->detectIndexPatterns();
        }

        inline Argument<Base> makeGraphClones(const EquationPattern<Base>& eq,
                                              OperationNode<Base>& node) {

            assert(node.getOperationType() != CGInvOp);

            size_t id = node.getVariableID();

            if (id > 0) {
                // been here before
                return Argument<Base>(*clonesTemporary_.at(id));
            }

            node.increaseUsageCount();

            id = idCounter_++;
            node.setVariableID(id);

            if (node.getColor() > 0 || node.getOperationType() == CGArrayCreationOp) {
                /**
                 * part of the operation path that depends on the loop indexes
                 * or its an array with constant elements
                 */
                const std::vector<Argument<Base> >& args = node.getArguments();
                size_t arg_size = args.size();
                std::vector<Argument<Base> > cloneArgs(arg_size);

                for (size_t a = 0; a < arg_size; a++) {
                    OperationNode<Base>* argOp = args[a].getOperation();
                    if (argOp == NULL) {
                        // parameter
                        cloneArgs[a] = Argument<Base>(*args[a].getParameter());
                    } else {
                        // variable
                        if (argOp->getOperationType() == CGInvOp) {
                            cloneArgs[a] = Argument<Base>(createIndependentClone(&node, a, *argOp));
                        } else {
                            cloneArgs[a] = makeGraphClones(eq, *argOp);
                        }
                    }
                }

                OperationNode<Base>* cloneOp = new OperationNode<Base>(
                        node.getOperationType(),
                        node.getInfo(),
                        cloneArgs);
                handler_.manageOperationNode(cloneOp);

                clonesTemporary_[id] = cloneOp;
                return Argument<Base>(*cloneOp);

            } else {
                /**
                 * temporary variable used in all iterations
                 * (does not depend on indexes)
                 */
                return Argument<Base>(makeTemporaryVarClone(node));
            }
        }

        inline OperationNode<Base>& createIndependentClone(OperationNode<Base>* operation,
                                                           size_t argumentIndex,
                                                           OperationNode<Base>& independent) {

            // is it an indexed independent?
            typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::const_iterator it = indexedOpIndep.op2Arguments.find(operation);
            if (it != indexedOpIndep.op2Arguments.end()) {
                const OperationIndexedIndependents<Base>& yyy = it->second;

                if (!yyy.arg2Independents[argumentIndex].empty()) {
                    // yes
                    return getIndexedIndependentClone(operation, argumentIndex);
                }
            }

            // it is constant for all operation
            return getNonIndexedIndependentClone(independent);
        }

        OperationNode<Base>& getIndexedIndependentClone(const OperationNode<Base>* operation,
                                                        size_t argIndex) {
            assert(operation == NULL || operation->getArguments().size() > argIndex);
            assert(operation == NULL || operation->getArguments()[argIndex].getOperation() != NULL);
            assert(operation == NULL || operation->getArguments()[argIndex].getOperation()->getOperationType() == CGInvOp);

            OperationArgumentsIndepOrder<Base>* args2Order = op2Arg2IndepOrder_.at(operation);
            IndependentOrder<Base>* indepOrder = args2Order->arg2Order.at(argIndex);

            typename std::map<const IndependentOrder<Base>*, OperationNode<Base>*>::const_iterator it;
            it = indexedIndep2clone_.find(indepOrder);
            if (it != indexedIndep2clone_.end()) {
                return *it->second;
            } else {
                CG<Base> newIndep;
                handler_.makeVariable(newIndep);
                independentsIndexed_[newIndep.getOperationNode()] = IndexValue(nIndependents_, newIndep);
                nIndependents_++;

                OperationNode<Base>* clone = newIndep.getOperationNode();
                indexedIndep2clone_[indepOrder] = clone;
                return *clone;
            }
        }

        OperationNode<Base>& getNonIndexedIndependentClone(const OperationNode<Base>& node) {
            assert(node.getOperationType() == CGInvOp);

            typename std::map<const OperationNode<Base>*, OperationNode<Base>*>::iterator it;
            it = orig2ConstIndepClone_.find(&node);
            if (it != orig2ConstIndepClone_.end()) {
                return *it->second;
            }

            CG<Base> newIndep;
            handler_.makeVariable(newIndep);
            independentsNonIndexed_[newIndep.getOperationNode()] = IndexValue(nIndependents_, newIndep);
            nIndependents_++;

            OperationNode<Base>* clone = newIndep.getOperationNode();
            orig2ConstIndepClone_[&node] = clone;
            return *clone;
        }

        /**
         * Creates a temporary variable that does NOT depend on the loop indexes
         * 
         * @param node The original node
         * @return the clone
         */
        OperationNode<Base>& makeTemporaryVarClone(OperationNode<Base>& node) {
            assert(node.getOperationType() != CGInvOp);
            assert(node.getOperationType() != CGArrayCreationOp);

            CG<Base> newIndep;
            handler_.makeVariable(newIndep);
            OperationNode<Base>* cloneOp = newIndep.getOperationNode();

            temporaryClone2Orig_[cloneOp] = &node;
            independentsTemp_[cloneOp] = IndexValue(nIndependents_, newIndep);
            nIndependents_++;

            size_t id = idCounter_++;
            clonesTemporary_[id] = cloneOp;
            node.setVariableID(id);

            return *cloneOp;
        }

        Loop(const Loop<Base>& other); // not implemented
        Loop& operator=(const Loop<Base>& rhs); // not implemented

        static void resetCounters(OperationNode<Base>& node) {
            if (node.getVariableID() == 0) {
                return;
            }

            node.setVariableID(0);
            node.resetUsageCount();
            node.setLastUsageEvaluationOrder(0);

            const std::vector<Argument<Base> >& args = node.getArguments();
            size_t arg_size = args.size();
            for (size_t i = 0; i < arg_size; i++) {
                if (args[i].getOperation() != NULL) {
                    resetCounters(*args[i].getOperation());
                }
            }
        }

        /**
         * structure used to sort the loop's indexed independent variables
         */
        struct IndexedIndepSorter {
            const std::map<const OperationNode<Base>*, const IndependentOrder<Base>*>& clone2indexedIndep;
            const std::map<const OperationNode<Base>*, size_t>& origModelIndepOrder;

            IndexedIndepSorter(const std::map<const OperationNode<Base>*, const IndependentOrder<Base>*>& clone2indexedIndep_,
                               const std::map<const OperationNode<Base>*, size_t>& origModelIndepOrder_) :
                clone2indexedIndep(clone2indexedIndep_),
                origModelIndepOrder(origModelIndepOrder_) {
            }

            bool operator()(const OperationNode<Base>* node1,
                    const OperationNode<Base>* node2) {
                const IndependentOrder<Base>* indepOrder1 = clone2indexedIndep.at(node1);
                const IndependentOrder<Base>* indepOrder2 = clone2indexedIndep.at(node2);
                assert(indepOrder1->order.size() == indepOrder2->order.size());

                size_t size = indepOrder1->order.size();
                for (size_t j = 0; j < size; j++) {
                    const OperationNode<Base>* indep1 = indepOrder1->order[j];
                    const OperationNode<Base>* indep2 = indepOrder2->order[j];
                    // some variables are not used in all iterations
                    if (indep1 == NULL) {
                        if (indep2 == NULL) {
                            continue;
                        }
                        return false;
                    } else if (indep2 == NULL) {
                        return true;
                    }

                    size_t index1 = origModelIndepOrder.at(indep1);
                    size_t index2 = origModelIndepOrder.at(indep2);
                    if (index1 < index2)
                        return true;
                    else if (index1 > index2)
                        return false;
                }

                assert(false); // should never get here
                return false;
            }
        };

    };

}

#endif
