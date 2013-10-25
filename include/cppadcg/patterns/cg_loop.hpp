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
         * The number of iterations this loop will have
         */
        const size_t iterationCount;
        /**
         * The equations inside the loop
         */
        std::set<EquationPattern<Base>*> equations;
        /**
         * Which argument positions of operations (from the reference dependent)
         * use indexed independent variables
         */
        IndexedIndependent<Base> indexedOpIndep;
        /**
         * Which dependents (equation indexes) must be evaluated at the same 
         * time due to shared temporary variables
         */
        std::vector<std::set<size_t> > linkedDependents;
    private:
        /**
         * code handler for the operations in the loop
         */
        LoopCodeHandler<Base> handler_;
        /**
         * The order of the equation patterns in the loop's tape
         */
        std::map<EquationPattern<Base>*, size_t> equationOrder_;
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
         * The evaluated dependents in each loop iteration
         * ([iteration 1]{dep1, dep3, ...}; [iteration 2]{dep5, dep6, ...}; ...)
         */
        std::vector<std::set<size_t> > iterationDependents_;
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

        Loop(EquationPattern<Base>& eq) :
            iterationCount(eq.dependents.size()),
            nIndependents_(0),
            loopModel_(NULL) {
            //indexedOpIndep(eq.indexedOpIndep), // must be determined later for a different reference
            //constOperationIndependents(eq.constOperationIndependents) {
            equations.insert(&eq);
        }

        void addEquation(EquationPattern<Base>& eq,
                         const std::vector<std::set<size_t> >& newLinkedDependents) {
            assert(eq.dependents.size() == iterationCount);

            if (equations.insert(&eq).second)
                return; //it had already been added

            //indexedOpIndep.combine(eq.indexedOpIndep); // must be determined later for a different reference
            //constOperationIndependents.insert(eq.constOperationIndependents.begin(),
            //                                  eq.constOperationIndependents.end());
        }

        const std::vector<std::set<size_t> >& getIterationDependents() const {
            return iterationDependents_;
        }

        LoopModel<Base>* getModel() const {
            return loopModel_;
        }

        LoopModel<Base>* releaseLoopModel() {
            LoopModel<Base>* loopAtomic = loopModel_;
            loopModel_ = NULL;
            return loopAtomic;
        }

        /**
         * Attempts to combine the provided loop with the current one
         * 
         * @param other The other loop
         * @return true if the other loop was merged into this one, false otherwise
         */
        bool merge(Loop<Base>& other) {
            if (iterationCount != other.iterationCount) {
                return false;
            }

            size_t osize = other.linkedDependents.size();
            for (size_t oe = 0; oe < osize; oe++) {
                combineDependentRelations(other.linkedDependents, linkedDependents);
            }

            equations.insert(other.equations.begin(), other.equations.end());
            other.equations.clear(); // so that it does not delete the equations

            return true;
        }

        /**
         * Checks if there would be only one depenedent from the same equation
         * for the same loop iteration
         * 
         * @param newRelations the indexes of the dependent variables that must
         *                     evaluated in the same iteration
         * @param dep2Eq maps the dependent variable indexes to their equation patterns
         * @return true if there is only dependent from each equation in all 
         *         loop iterations, false otherwise
         */
        bool isCompatibleWithDepRelationships(const std::vector<std::set<size_t> >& newRelations,
                                              const std::map<size_t, EquationPattern<Base>*>& dep2Eq) const {
            /**
             * create what would be the new relationships
             */
            std::vector<std::set<size_t> > linksAfter = linkedDependents; //copy

            size_t osize = newRelations.size();
            for (size_t oe = 0; oe < osize; oe++) {
                combineDependentRelations(newRelations[oe], linksAfter);
            }

            /**
             * Cannot have more than one dependent from the same equation in
             * the same relationship
             */
            size_t l_size = linksAfter.size();
            for (size_t l = 0; l < l_size; l++) {
                const std::set<size_t>& deps = linksAfter[l];
                std::set<EquationPattern<Base>*> equations;

                std::set<size_t>::const_iterator it;
                for (it = deps.begin(); it != deps.end(); ++it) {
                    EquationPattern<Base>* eq = dep2Eq.at(*it);
                    if (!equations.insert(eq).second) {
                        return false; // more than one dependent for the same equation
                    }
                }
            }

            return true; // all OK
        }

        void createLoopModel(const std::vector<CG<Base> >& dependents,
                             const std::vector<CG<Base> >& independents,
                             const std::map<size_t, EquationPattern<Base>*>& dep2Equation,
                             std::map<OperationNode<Base>*, size_t>& origTemp2Index) throw (CGException) {

            generateDependentLoopIndexes(dep2Equation);

            const std::set<size_t>& firstItDep = iterationDependents_[0];
            assert(firstItDep.size() == equations.size());

            std::set<size_t>::const_iterator itDep;
            for (itDep = firstItDep.begin(); itDep != firstItDep.end(); ++itDep) {
                size_t dep = *itDep;
                EquationPattern<Base>::uncolor(dependents[dep].getOperationNode());
            }

            for (itDep = firstItDep.begin(); itDep != firstItDep.end(); ++itDep) {
                size_t dep = *itDep;
                EquationPattern<Base>* eq = dep2Equation.at(dep);

                // operations that use indexed independent variables in the first iteration
                std::set<const OperationNode<Base>*> indexedOperations;

                eq->colorIndexedPath(dep, dependents, 1, indexedOperations);
                if (dep == eq->depRefIndex) {
                    indexedOpIndep.op2Arguments.insert(eq->indexedOpIndep.op2Arguments.begin(),
                                                       eq->indexedOpIndep.op2Arguments.end());
                } else {
                    generateLoopOperationReferences(eq, indexedOperations);
                }
            }

            generateIndependentLoopIndexes(dep2Equation);

            resetCounters(dependents);

            createLoopTapeNModel(dependents, independents, dep2Equation, origTemp2Index);

            /**
             * Clean-up
             */
            for (itDep = firstItDep.begin(); itDep != firstItDep.end(); ++itDep) {
                size_t dep = *itDep;
                EquationPattern<Base>::uncolor(dependents[dep].getOperationNode());
            }
            resetCounters(dependents);
        }

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
         *                        public static methods
         **********************************************************************/
        static bool canCombineEquations(EquationPattern<Base>& eq1,
                                        EquationPattern<Base>& eq2,
                                        const OperationNode<Base>& sharedTemp) {
            // convert to the reference operation
            OperationNode<Base>* sharedTempRef1 = eq1.operationEO2Reference.at(sharedTemp.getEvaluationOrder());

            // must have indexed independents at the same locations in all equations
            const std::set<const OperationNode<Base>*> opWithIndepArgs = EquationPattern<Base>::findOperationsUsingIndependents(*sharedTempRef1);

            // must have indexed independents at the same locations in both equations
            typename std::set<const OperationNode<Base>*>::const_iterator itOp;
            for (itOp = opWithIndepArgs.begin(); itOp != opWithIndepArgs.end(); ++itOp) {
                const OperationNode<Base>* op1 = *itOp;

                typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base>*>::const_iterator indexed1It;
                indexed1It = eq1.indexedOpIndep.op2Arguments.find(op1);

                typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base>*>::const_iterator indexed2It;
                indexed2It = eq2.indexedOpIndep.op2Arguments.find(eq2.operationEO2Reference.at(op1->getEvaluationOrder()));

                if (indexed1It == eq1.indexedOpIndep.op2Arguments.end()) {
                    if (indexed2It != eq2.indexedOpIndep.op2Arguments.end()) {
                        return false;
                    }
                } else {
                    if (indexed2It == eq2.indexedOpIndep.op2Arguments.end()) {
                        return false;
                    }
                    const OperationIndexedIndependents<Base>& indexed1Ops = *indexed1It->second;
                    const OperationIndexedIndependents<Base>& indexed2Ops = *indexed2It->second;
                    if (indexed1Ops.arg2Independents.size() != indexed2Ops.arg2Independents.size()) {
                        return false;
                    }

                    typename std::map<size_t, std::map<size_t, const OperationNode<Base>*> >::const_iterator argPos1It = indexed1Ops.arg2Independents.begin();
                    typename std::map<size_t, std::map<size_t, const OperationNode<Base>*> >::const_iterator argPos2It = indexed2Ops.arg2Independents.begin();

                    for (; argPos2It != indexed2Ops.arg2Independents.end(); ++argPos1It, ++argPos2It) {
                        if (argPos1It->first != argPos2It->first)
                            return false;
                    }
                }
            }

            return true;
        }

        static void combineDependentRelations(const std::vector<std::set<size_t> >& newDependentLinks,
                                              std::vector<std::set<size_t> >& linkedDependents) {
            size_t size = newDependentLinks.size();
            for (size_t l = 0; l < size; l++) {
                combineDependentRelations(newDependentLinks[l], linkedDependents);
            }
        }

        static void combineDependentRelations(const std::set<size_t>& oSameIndex,
                                              std::vector<std::set<size_t> >& linkedDependents) {
            std::set<size_t> pos;

            if (linkedDependents.size() > 0) {
                std::set<size_t>::const_iterator oIt;
                for (oIt = oSameIndex.begin(); oIt != oSameIndex.end(); ++oIt) {
                    size_t otherDep = *oIt;
                    long e = findLinkedDependent(otherDep, linkedDependents);
                    if (e >= 0) {
                        pos.insert(e); // found
                    }
                }
            }

            if (pos.size() > 0) {
                size_t init = *pos.begin();
                linkedDependents[init].insert(oSameIndex.begin(), oSameIndex.end());

                std::set<size_t>::const_reverse_iterator itr;
                for (itr = pos.rbegin(); itr != pos.rend();) {
                    size_t e = *itr;
                    ++itr;
                    if (e != init) {
                        linkedDependents[init].insert(linkedDependents[e].begin(), linkedDependents[e].end());
                        linkedDependents.erase(linkedDependents.begin() + e);
                    }
                }

            } else {
                size_t size = linkedDependents.size();
                linkedDependents.resize(size + 1);
                linkedDependents[size].insert(oSameIndex.begin(), oSameIndex.end());
            }
        }

        /***********************************************************************
         *                            private
         **********************************************************************/
    private:

        inline long findLinkedDependent(size_t dep) const {
            return findLinkedDependent(dep, linkedDependents);
        }

        static inline long findLinkedDependent(size_t dep,
                                               const std::vector<std::set<size_t> >& linkedDependents) {
            size_t size = linkedDependents.size();
            for (size_t pos = 0; pos < size; pos++) {
                const std::set<size_t>& sameIndex = linkedDependents[pos];
                if (sameIndex.find(dep) != sameIndex.end()) {
                    return pos;
                }
            }

            return -1;
        }

        inline bool containsDepOutsideIndex(size_t dep, size_t el) const {
            return containsDepOutsideIndex(dep, el, linkedDependents);
        }

        static inline bool containsDepOutsideIndex(size_t dep,
                                                   size_t pos,
                                                   const std::vector<std::set<size_t> >& linkedDependents) {
            size_t size = linkedDependents.size();
            for (size_t e = 0; e < size; e++) {
                if (e != pos) {
                    const std::set<size_t>& sameIndex = linkedDependents[e];
                    if (sameIndex.find(dep) != sameIndex.end()) {
                        // already used with a different dependent (index)
                        return true;
                    }
                }
            }

            return false;
        }

        void generateDependentLoopIndexes(const std::map<size_t, EquationPattern<Base>*>& dep2Equation) {
            iterationDependents_.clear();
            iterationDependents_.resize(iterationCount);

            equationOrder_.clear();

            /**
             * assign a dependent variable from each equation to an iteration
             */
            std::map<EquationPattern<Base>*, std::set<size_t> > depsInEq;

            typename std::set<EquationPattern<Base>*>::const_iterator eqIt;
            for (eqIt = equations.begin(); eqIt != equations.end(); ++eqIt) {
                EquationPattern<Base>* eq = *eqIt;
                depsInEq[eq] = eq->dependents;
            }

            /**
             * Determine the dependents of the first iteration
             */
            std::set<EquationPattern<Base>* > eqs = equations;
            std::set<size_t>& iterationDeps0 = iterationDependents_[0];
            std::map<EquationPattern<Base>*, std::vector<size_t> > masterEquations; // equation that define the order

            while (eqs.size() > 0) {
                // find the lowest dependent variable index of the remaining equations
                size_t lowestDep = std::numeric_limits<size_t>::max();
                typename std::set<EquationPattern<Base>* >::iterator itEq;
                for (itEq = eqs.begin(); itEq != eqs.end(); ++itEq) {
                    EquationPattern<Base>* eq = *itEq;
                    size_t lDep = *depsInEq.at(eq).begin();
                    if (lDep < lowestDep)
                        lowestDep = lDep;
                }

                EquationPattern<Base>* masterEq = dep2Equation.at(lowestDep);
                std::vector<size_t>& depOrder = masterEquations[masterEq];
                depOrder.insert(depOrder.begin(), masterEq->dependents.begin(), masterEq->dependents.end());

                long pos = findLinkedDependent(lowestDep);
                if (pos >= 0) {
                    // assign the dependent to the first iteration with all its relationships
                    std::set<size_t>::const_iterator itDep;
                    for (itDep = linkedDependents[pos].begin(); itDep != linkedDependents[pos].end(); ++itDep) {
                        size_t dep = *itDep;
                        iterationDeps0.insert(dep);
                        EquationPattern<Base>* eq = dep2Equation.at(dep);
                        size_t eqo_size = equationOrder_.size();
                        equationOrder_[eq] = eqo_size;
                        depsInEq.at(eq).erase(dep);
                        eqs.erase(eq);
                    }
                } else {
                    iterationDeps0.insert(lowestDep);
                    size_t eqo_size = equationOrder_.size();
                    equationOrder_[masterEq] = eqo_size;
                    depsInEq.at(masterEq).erase(lowestDep);
                    eqs.erase(masterEq);
                }
            }

            /**
             * Determine the dependent variable indexes of the following 
             * iterations
             */
            for (size_t i = 1; i < iterationCount; i++) {
                std::set<size_t>& iterationDeps = iterationDependents_[i];

                typename std::map<EquationPattern<Base>*, std::vector<size_t> >::const_iterator itm;
                for (itm = masterEquations.begin(); itm != masterEquations.end(); ++itm) {
                    const std::vector<size_t>& depOrder = itm->second;
                    size_t mDep = depOrder[i];
                    long pos = findLinkedDependent(mDep);
                    if (pos >= 0) {
                        // place the dependent with all its relationships
                        std::set<size_t>::const_iterator itDep;
                        for (itDep = linkedDependents[pos].begin(); itDep != linkedDependents[pos].end(); ++itDep) {
                            iterationDeps.insert(*itDep);
                        }
                    } else {
                        iterationDeps.insert(mDep);
                    }
                }
            }

        }

        void generateLoopOperationReferences(const EquationPattern<Base>* eq,
                                             const std::set<const OperationNode<Base>*>& indexedOperations) {
            typename std::set<const OperationNode<Base>*>::const_iterator it;
            for (it = indexedOperations.begin(); it != indexedOperations.end(); ++it) {
                const OperationNode<Base>* opLoopRef = *it;
                const OperationNode<Base>* opEqRef = eq->operationEO2Reference.at(opLoopRef->getEvaluationOrder());
                indexedOpIndep.op2Arguments[opLoopRef] = eq->indexedOpIndep.arguments(opEqRef);
            }
        }

        void generateIndependentLoopIndexes(const std::map<size_t, EquationPattern<Base>*>& dep2Equation) {
            /**
             * find indexed independents with the same order
             * so that they can be defined as a single variable 
             * in the atomic function
             */
            std::map<size_t, size_t> dep2Iteration;
            for (size_t iter = 0; iter < iterationCount; iter++) {
                const std::set<size_t>& deps = iterationDependents_[iter];

                std::set<size_t>::const_iterator itDeps;
                for (itDeps = deps.begin(); itDeps != deps.end(); ++itDeps) {
                    dep2Iteration[*itDeps] = iter;
                }
            }

            // loop all operations from the reference dependents which use indexed independents
            typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base>*>::const_iterator it;
            for (it = indexedOpIndep.op2Arguments.begin(); it != indexedOpIndep.op2Arguments.end(); ++it) {
                const OperationNode<Base>* operation = it->first;
                const OperationIndexedIndependents<Base>& opInd = *it->second;

                OperationArgumentsIndepOrder<Base>* arg2orderPos = new OperationArgumentsIndepOrder<Base>();
                op2Arg2IndepOrder_[operation] = arg2orderPos;

                arg2IndepOrder_.push_front(arg2orderPos);

                // loop all arguments
                typename std::map<size_t, std::map<size_t, const OperationNode<Base>* > >::const_iterator it2;
                for (it2 = opInd.arg2Independents.begin(); it2 != opInd.arg2Independents.end(); ++it2) {
                    size_t argumentIndex = it2->first;
                    const std::map<size_t, const OperationNode<Base>*>& dep2Indep = it2->second;

                    std::vector<const OperationNode<Base>*> order(iterationCount);

                    /**
                     * create the independent variable order
                     */
                    assert(dep2Indep.size() == iterationCount);
                    typename std::map<size_t, const OperationNode<Base>*>::const_iterator itDep2Indep;
                    for (itDep2Indep = dep2Indep.begin(); itDep2Indep != dep2Indep.end(); ++itDep2Indep) {
                        size_t dep = itDep2Indep->first;
                        const OperationNode<Base>* indep = itDep2Indep->second;

                        // find one of the dependents -> iteration index
                        size_t iterationIndex = dep2Iteration.at(dep);

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
                        for (size_t iterationIndex = 0; iterationIndex < iterationCount; iterationIndex++) {
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
         * Resets the counters for the first loop iteration
         * 
         * @param dependents
         */
        void resetCounters(const std::vector<CG<Base> >& dependents) {
            std::set<size_t>::const_iterator depIt;
            for (depIt = iterationDependents_[0].begin(); depIt != iterationDependents_[0].end(); ++depIt) {
                size_t depIndex = *depIt;
                OperationNode<Base>* node = dependents[depIndex].getOperationNode();
                if (node != NULL) {
                    Loop<Base>::resetCounters(*node);
                }
            }
        }

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
             * create the new/clone operations for the first iteration only
             */
            nIndependents_ = 0;
            idCounter_ = 1;

            std::vector<CGB> deps(iterationDependents_[0].size());

            std::set<size_t>::const_iterator depIt;
            for (depIt = iterationDependents_[0].begin(); depIt != iterationDependents_[0].end(); ++depIt) {
                size_t depIndex = *depIt;
                EquationPattern<Base>* eq = dep2Equation.at(depIndex);
                OperationNode<Base>* node = dependents[depIndex].getOperationNode();

                Argument<Base> aClone;

                if (node != NULL) {
                    if (node->getOperationType() == CGInvOp) {
                        OperationNode<Base>* nodeRef = NULL;
                        aClone = Argument<Base>(createIndependentClone(*eq, NULL, 0, *node, nodeRef));
                    } else {
                        aClone = makeGraphClones(*eq, *node);
                    }
                } else {
                    aClone = Argument<Base>(dependents[depIndex].getValue());
                }

                size_t i = equationOrder_.at(eq);
                deps[i] = CGB(handler_, aClone);
            }

            /*******************************************************************
             * create the tape for the first iteration
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
             * Reorder independent variables for the new tape (1st iteration only)
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
            std::vector<std::vector<size_t> > dependentOrigIndexes(equations.size(), std::vector<size_t> (iterationCount));
            for (size_t it = 0; it < iterationCount; it++) {
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
            std::vector<std::vector<size_t> > indexedIndependents(nIndexed, std::vector<size_t>(iterationCount));
            for (size_t j = 0; j < indexedCloneOrder.size(); j++) {
                const IndependentOrder<Base>* origOrder = clone2indexedIndep.at(indexedCloneOrder[j]);
                for (size_t it = 0; it < iterationCount; it++) {
                    size_t index = origModelIndepOrder.at(origOrder->order[it]);
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
                    iterationCount,
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
                OperationNode<Base>* nodeRef = NULL;
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
                            cloneArgs[a] = Argument<Base>(createIndependentClone(eq, &node, a, *argOp, nodeRef));
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

        inline OperationNode<Base>& createIndependentClone(const EquationPattern<Base>& eq,
                                                           OperationNode<Base>* operation, size_t argumentIndex,
                                                           OperationNode<Base>& independent,
                                                           OperationNode<Base>*& operationRef) {
            if (operationRef == NULL && operation != NULL) {
                operationRef = eq.operationEO2Reference.at(operation->getEvaluationOrder());
            }

            // is it an indexed independent?
            OperationIndexedIndependents<Base>* yyy = eq.indexedOpIndep.find(operationRef);
            if (yyy != NULL) {
                typename std::map<size_t, std::map<size_t, const OperationNode<Base>*> >::const_iterator itadi;
                itadi = yyy->arg2Independents.find(argumentIndex);
                if (itadi != yyy->arg2Independents.end()) {
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
                    size_t index1 = origModelIndepOrder.at(indepOrder1->order[j]);
                    size_t index2 = origModelIndepOrder.at(indepOrder2->order[j]);
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
