#ifndef CPPAD_CG_DEPENDENT_PATTERN_MATCHER_INCLUDED
#define CPPAD_CG_DEPENDENT_PATTERN_MATCHER_INCLUDED
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

//#define CPPADCG_PRINT_DEBUG

namespace CppAD {

    template<class Base>
    class UniqueEquationPair {
    public:
        EquationPattern<Base>* eq1;
        EquationPattern<Base>* eq2;
    public:

        inline UniqueEquationPair(EquationPattern<Base>* e1, EquationPattern<Base>* e2) {
            if (e1->depRefIndex < e2->depRefIndex) {
                eq1 = e1;
                eq2 = e2;
            } else {
                eq1 = e2;
                eq2 = e1;
            }
        }
    };

    template<class Base>
    inline bool operator<(const UniqueEquationPair<Base>& p1, const UniqueEquationPair<Base>& p2) {
        return p1.eq1->depRefIndex < p2.eq1->depRefIndex || (!(p2.eq1->depRefIndex < p1.eq1->depRefIndex) && p1.eq2->depRefIndex < p2.eq2->depRefIndex);
    }

    /**
     * Finds common patterns in operation graphs
     */
    template<class Base>
    class DependentPatternMatcher {
    public:
        typedef CG<Base> CGBase;
    private:

        enum INDEXED_OPERATION_TYPE {
            INDEXED_OPERATION_TYPE_INDEXED,
            INDEXED_OPERATION_TYPE_NONINDEXED,
            INDEXED_OPERATION_TYPE_BOTH
        };
        typedef std::pair<INDEXED_OPERATION_TYPE, size_t> Indexed2OpCountType;
        typedef std::map<size_t, std::map<size_t, std::map<OperationNode<Base>*, Indexed2OpCountType> > > Dep1Dep2SharedType;
        typedef std::pair<size_t, size_t> DepPairType;
        typedef std::map<size_t, std::map<DepPairType, const std::map<OperationNode<Base>*, Indexed2OpCountType>* > > TotalOps2validDepsType;
        typedef std::map<UniqueEquationPair<Base>, TotalOps2validDepsType*> Eq2totalOps2validDepsType;
        typedef std::map<size_t, Eq2totalOps2validDepsType> MaxOps2eq2totalOps2validDepsType;

    private:
        const std::vector<std::set<size_t> >& relatedDepCandidates_;
        std::vector<CGBase> dependents_; // a copy
        std::vector<CGBase>& independents_;
        std::vector<EquationPattern<Base>*> equations_;
        EquationPattern<Base>* eqCurr_;
        std::map<size_t, EquationPattern<Base>*> dep2Equation_;
        std::map<EquationPattern<Base>*, Loop<Base>*> equation2Loop_;
        std::vector<Loop<Base>*> loops_;
        /**
         * Equations which cannot be in the same loop
         */
        std::map<EquationPattern<Base>*, std::set<EquationPattern<Base>*> > incompatible_;
        /**
         * 
         */
        std::map<UniqueEquationPair<Base>, Dep1Dep2SharedType> equationShared_;
        /**
         * maps the original model nodes used as temporary non-indexed variables
         * by the loops to an index k
         */
        std::map<OperationNode<Base>*, size_t> origTemp2Index_;
        std::vector<std::set<size_t> > id2Deps;
        size_t idCounter_;
        /// used to mark visited nodes and indexed nodes
        size_t color_;
    public:

        /**
         * Creates a new DependentPatternMatcher
         * 
         * @param relatedDepCandidates Groups of dependent variable indexes that
         *                             are believed to have the same expression
         *                             pattern.
         * @param dependents The dependent variable values
         * @param independents The independent variable values
         */
        DependentPatternMatcher(const std::vector<std::set<size_t> >& relatedDepCandidates,
                                const std::vector<CGBase>& dependents,
                                std::vector<CGBase>& independents) :
            relatedDepCandidates_(relatedDepCandidates),
            dependents_(dependents),
            independents_(independents),
            idCounter_(0),
            color_(0) {
            CPPADCG_ASSERT_UNKNOWN(independents_.size() > 0);
            CPPADCG_ASSERT_UNKNOWN(independents_[0].getCodeHandler() != NULL);
            equations_.reserve(relatedDepCandidates_.size());
        }

        const std::vector<EquationPattern<Base>*>& getEquationPatterns() const {
            return equations_;
        }

        const std::vector<Loop<Base>*>& getLoops() const {
            return loops_;
        }

        /**
         * Detects common equation patterns and generates a new tape for the
         * model using loops.
         * This method should only be called once!
         * 
         * @param nonLoopTape The new tape without the loops or NULL if there
         *                    are no non-indexed expressions in the model
         * @param loopTapes The models for each loop (must be deleted by the user)
         */
        virtual void generateTapes(LoopFreeModel<Base>*& nonLoopTape,
                                   std::set<LoopModel<Base>*>& loopTapes) throw (CGException) {

            for (size_t j = 0; j < independents_.size(); j++) {
                std::vector<size_t>& info = independents_[j].getOperationNode()->getInfo();
                info.resize(1);
                info[0] = j;
            }

            findLoops();

            nonLoopTape = createNewTape();

            loopTapes.clear();
            for (size_t l = 0; l < loops_.size(); l++) {
                Loop<Base>* loop = loops_[l];
                loopTapes.insert(loop->releaseLoopModel());
            }
        }

        virtual ~DependentPatternMatcher() {
            for (size_t l = 0; l < loops_.size(); l++) {
                delete loops_[l];
            }
        }

    private:

        /**
         * Attempts to detect common patterns in the equations and generate 
         * loops from these patterns.
         * 
         * @return information about the detected loops
         */
        virtual std::vector<Loop<Base>*> findLoops() throw (CGException) {
            using namespace std;

            size_t rSize = relatedDepCandidates_.size();
            for (size_t r = 0; r < rSize; r++) {
                const std::set<size_t>& candidates = relatedDepCandidates_[r];
                std::set<size_t>::const_iterator itDep;
                for (itDep = candidates.begin(); itDep != candidates.end(); ++itDep) {
                    size_t iDep = *itDep;
                    OperationNode<Base>* node = dependents_[iDep].getOperationNode();
                    if (node != NULL && node->getOperationType() == CGInvOp) {
                        /**
                         * indexed/nonindexed independents are marked at the 
                         * operation that uses them, therefore currently there 
                         * is no way to make a distinction between
                         *     yi = xi and y_(i+1) = x_(i+1)
                         * since both operations which use indexed independents
                         * have no operation, consequently an alias is used
                         */
                        CodeHandler<Base>* handler = dependents_[iDep].getCodeHandler();
                        dependents_[iDep] = handler->createCG(new OperationNode<Base>(CGAliasOp, Argument<Base>(*node)));
                    }
                }
            }

            // assign a unique Id to each node
            assignIds();
            id2Deps.resize(idCounter_ + 1);

            /**
             * Determine the equation patterns
             */
            findRelatedVariables();

            const size_t eq_size = equations_.size();
            for (size_t e = 0; e < eq_size; e++) {
                EquationPattern<Base>* eq = equations_[e];

                set<size_t>::const_iterator depIt;
                for (depIt = eq->dependents.begin(); depIt != eq->dependents.end(); ++depIt) {
                    dep2Equation_[*depIt] = eq;
                }
            }

            loops_.reserve(eq_size);

            SmartSetPointer<set<size_t> > dependentRelations;
            std::vector<set<size_t>*> dep2Relations(dependents_.size(), NULL);
            map<size_t, set<size_t> > dependentBlackListRelations;

            /*******************************************************************
             * Combine related equations in the same loops
             * (equations that share temporary variables)
             ******************************************************************/
            /**
             * Find and organize relationships
             */
            color_++; // this is the color used to mark indexed nodes (must be higher than any previously used color)

            for (size_t e = 0; e < eq_size; e++) {
                EquationPattern<Base>* eq = equations_[e];
                eqCurr_ = eq;
                set<size_t>::const_iterator depIt;

                for (depIt = eq->dependents.begin(); depIt != eq->dependents.end(); ++depIt) {
                    OperationNode<Base>* node = dependents_[*depIt].getOperationNode();
                    // will define the dependents associated with each operation
                    markOperationsWithDependent(node, *depIt);
                }

                /**
                 * Find shared operations with the previous equation patterns
                 */
                if (e > 0) {
                    for (depIt = eq->dependents.begin(); depIt != eq->dependents.end(); ++depIt) {
                        findSharedTemporaries(dependents_[*depIt], *depIt); // a color is used to mark indexed paths
                    }

                    /**
                     * clean-up
                     */
                    for (depIt = eq->dependents.begin(); depIt != eq->dependents.end(); ++depIt) {
                        OperationNode<Base>* node = dependents_[*depIt].getOperationNode();
                        EquationPattern<Base>::uncolor(node); // must uncolor
                        EquationPattern<Base>::clearUsageCount(node); // must reset usage count
                    }
                }

                // create a loop for this equation
                Loop<Base>* loop = new Loop<Base>(*eq);
                loops_.push_back(loop);
                equation2Loop_[eq] = loop;
            }

            /*******************************************************************
             * Attempt to combine loops with shared variables
             ******************************************************************/
            MaxOps2eq2totalOps2validDepsType maxOps2Eq2totalOps2validDeps;
            Eq2totalOps2validDepsType eq2totalOps2validDeps;
            SmartListPointer<TotalOps2validDepsType> totalOps2validDepsMem;

            /**
             * First organize pairs of equation patterns with shared variables
             * by the maximum number of shared operations among two dependent
             * variables.
             */
            for (size_t l1 = 0; l1 < loops_.size(); l1++) {
                Loop<Base>* loop1 = loops_[l1];
                CPPADCG_ASSERT_UNKNOWN(loop1->equations.size() == 1);
                EquationPattern<Base>* eq1 = *loop1->equations.begin();

                for (size_t l2 = l1 + 1; l2 < loops_.size(); l2++) {
                    Loop<Base>* loop2 = loops_[l2];
                    CPPADCG_ASSERT_UNKNOWN(loop2->equations.size() == 1);
                    EquationPattern<Base>* eq2 = *loop2->equations.begin();

                    UniqueEquationPair<Base> eqRel(eq1, eq2);
                    typename map<UniqueEquationPair<Base>, Dep1Dep2SharedType>::const_iterator eqSharedit = equationShared_.find(eqRel);
                    if (eqSharedit == equationShared_.end())
                        continue; // nothing is shared between eq1 and eq2

                    const Dep1Dep2SharedType& dep1Dep2Shared = eqSharedit->second;

                    /**
                     * There are shared variables among the two equation patterns
                     */
                    TotalOps2validDepsType* totalOps2validDeps = new TotalOps2validDepsType();
                    totalOps2validDepsMem.push_back(totalOps2validDeps);
                    size_t maxOps = 0; // the maximum number of shared operations between two dependents

                    bool canCombine = true;

                    /***************************************************
                     * organize relations between dependents
                     **************************************************/
                    typename Dep1Dep2SharedType::const_iterator itDep1Dep2;
                    for (itDep1Dep2 = dep1Dep2Shared.begin(); itDep1Dep2 != dep1Dep2Shared.end(); ++itDep1Dep2) {
                        size_t dep1 = itDep1Dep2->first;
                        const map<size_t, map<OperationNode<Base>*, Indexed2OpCountType> >& dep2Shared = itDep1Dep2->second;

                        // multiple deps2 means multiple choices for a relation (only one dep1<->dep2 can be chosen)
                        typename map<size_t, map<OperationNode<Base>*, Indexed2OpCountType> >::const_iterator itDep2;
                        for (itDep2 = dep2Shared.begin(); itDep2 != dep2Shared.end(); ++itDep2) {
                            size_t dep2 = itDep2->first;
                            const map<OperationNode<Base>*, Indexed2OpCountType>& sharedTmps = itDep2->second;

                            size_t totalOps = 0; // the total number of operations performed by shared variables with dep2
                            typename map<OperationNode<Base>*, Indexed2OpCountType>::const_iterator itShared;
                            for (itShared = sharedTmps.begin(); itShared != sharedTmps.end(); ++itShared) {
                                if (itShared->second.first == INDEXED_OPERATION_TYPE_BOTH) {
                                    /**
                                     * one equation uses this temporary shared 
                                     * variable as an indexed variable while the 
                                     * other equation does not
                                     */
                                    canCombine = false;
                                    break;
                                } else {
                                    totalOps += itShared->second.second;
                                }
                            }

                            if (!canCombine) break;

                            DepPairType depRel(dep1, dep2);
                            (*totalOps2validDeps)[totalOps][depRel] = &sharedTmps;
                            maxOps = std::max(maxOps, totalOps);
                        }

                        if (!canCombine) break;
                    }

                    if (canCombine) {
                        maxOps2Eq2totalOps2validDeps[maxOps][eqRel] = totalOps2validDeps;
                        eq2totalOps2validDeps[eqRel] = totalOps2validDeps;
                    } else {
                        incompatible_[eq1].insert(eq2);
                        incompatible_[eq2].insert(eq1);
                        totalOps2validDepsMem.pop_back();
                        delete totalOps2validDeps;
                    }
                }
            }

            /**
             * Try to merge loops with shared variables 
             */
            typename MaxOps2eq2totalOps2validDepsType::const_reverse_iterator itMaxOps;
            for (itMaxOps = maxOps2Eq2totalOps2validDeps.rbegin(); itMaxOps != maxOps2Eq2totalOps2validDeps.rend(); ++itMaxOps) {
#ifdef CPPADCG_PRINT_DEBUG
                std::cout << "\n\nmaxOps: " << itMaxOps->first << "  count:" << itMaxOps->second.size() << std::endl;
#endif

                typename map<UniqueEquationPair<Base>, TotalOps2validDepsType*>::const_iterator itEqPair;
                for (itEqPair = itMaxOps->second.begin(); itEqPair != itMaxOps->second.end(); ++itEqPair) {
                    const UniqueEquationPair<Base>& eqRel = itEqPair->first;
#ifdef CPPADCG_PRINT_DEBUG
                    std::cout << "  eq1: " << *eqRel.eq1->dependents.begin() << "  eq2: " << *eqRel.eq2->dependents.begin() << std::endl;
#endif

                    Loop<Base>* loop1 = equation2Loop_.at(eqRel.eq1);
                    Loop<Base>* loop2 = equation2Loop_.at(eqRel.eq2);

                    if (loop1 == loop2)
                        continue; // already done
                    if (contains(incompatible_, eqRel.eq1, eqRel.eq2))
                        continue; // incompatible

                    /**
                     * backup so that it is possible to revert the new relations if
                     * required
                     */
                    SmartSetPointer<set<size_t> > dependentRelationsBak;
                    set<set<size_t>*>::const_iterator its;
                    for (its = dependentRelations.begin(); its != dependentRelations.end(); ++its) {
                        dependentRelationsBak.insert(new set<size_t>(**its));
                    }

                    // relationships between dependents for the resulting merged loop
                    set<set<size_t>*> loopRelations;

                    set<EquationPattern<Base>*> indexedLoopRelations;
                    std::vector<std::pair<EquationPattern<Base>*, EquationPattern<Base>*> > nonIndexedLoopRelations;

                    /**
                     * All equations from both loops must be compatible
                     */
                    bool compatible = isCompatible(loop1, loop2,
                                                   eq2totalOps2validDeps,
                                                   dep2Relations, dependentBlackListRelations, dependentRelations,
                                                   loopRelations, indexedLoopRelations, nonIndexedLoopRelations);

                    if (compatible) {
                        // merge the two loops

                        // update the loop of the equations
                        typename set<EquationPattern<Base>*>::const_iterator itle;
                        for (itle = loop2->equations.begin(); itle != loop2->equations.end(); ++itle) {
                            equation2Loop_[*itle] = loop1;
                        }
                        loop1->merge(*loop2, indexedLoopRelations, nonIndexedLoopRelations);

                        typename std::vector<Loop<Base>*>::iterator it = std::find(loops_.begin(), loops_.end(), loop2);
                        CPPADCG_ASSERT_UNKNOWN(it != loops_.end());
                        loops_.erase(it);
                        delete loop2;

                        loop1->setLinkedDependents(loopRelations);

                        // relation between loop1 and loop2 done!
                    } else {
                        // restore dependent relations
                        dependentRelations.s.swap(dependentRelationsBak.s);
                        // map each dependent to the relation set where it is present
                        std::fill(dep2Relations.begin(), dep2Relations.end(), (set<size_t>*) NULL);
                        for (its = dependentRelations.begin(); its != dependentRelations.end(); ++its) {
                            set<size_t>* relation = *its;
                            set<size_t>::const_iterator itd;
                            for (itd = relation->begin(); itd != relation->end(); ++itd) {
                                dep2Relations[*itd] = relation;
                            }
                        }

                    }

                }
            }

            /**
             * Determine the number of iterations in each loop
             */
            for (size_t l = 0; l < loops_.size(); l++) {
                loops_[l]->generateDependentLoopIndexes(dep2Equation_);
            }

            /*******************************************************************
             * Attempt to combine unrelated loops
             ******************************************************************/
            if (!loops_.empty()) {
                for (size_t l1 = 0; l1 < loops_.size() - 1; l1++) {
                    Loop<Base>* loop1 = loops_[l1];
                    for (size_t l2 = l1 + 1; l2 < loops_.size();) {
                        Loop<Base>* loop2 = loops_[l2];

                        bool canMerge = loop1->getIterationCount() == loop2->getIterationCount();
                        if (canMerge) {
                            // check if there are equations in the blacklist
                            canMerge = !find(loop1, loop2, incompatible_);
                        }

                        if (canMerge) {
                            loop1->mergeEqGroups(*loop2);
                            loops_.erase(loops_.begin() + l2);
                            delete loop2;
                        } else {
                            l2++;
                        }
                    }
                }
            }

            size_t l_size = loops_.size();

            /**
             * assign indexes (k) to temporary variables (non-indexed) used by loops
             */
            for (size_t l = 0; l < l_size; l++) {
                Loop<Base>* loop = loops_[l];

                //Generate a local model for the loop
                loop->createLoopModel(dependents_, independents_, dep2Equation_, origTemp2Index_);
            }

            /**
             * clean-up evaluation order
             */
            resetHandlerCounters();

            /**
             * clean-up colors
             */
            size_t m = dependents_.size();
            for (size_t i = 0; i < m; i++) {
                OperationNode<Base>* node = dependents_[i].getOperationNode();
                EquationPattern<Base>::uncolor(node);
            }

            return loops_;
        }

        /**
         * Determines whether or not two loops which shared temporary variables 
         * can be merged.
         * 
         * @param loop1 First loop
         * @param loop2 Second loop
         * @return true if the loops should be merged into a single loop
         */
        inline bool isCompatible(Loop<Base>* loop1,
                                 Loop<Base>* loop2,
                                 const Eq2totalOps2validDepsType& eq2totalOps2validDeps,
                                 std::vector<std::set<size_t>* >& dep2Relations,
                                 std::map<size_t, std::set<size_t> >& dependentBlackListRelations,
                                 SmartSetPointer<std::set<size_t> >& dependentRelations,
                                 std::set<std::set<size_t>*>& loopRelations,
                                 std::set<EquationPattern<Base>*>& indexedLoopRelations,
                                 std::vector<std::pair<EquationPattern<Base>*, EquationPattern<Base>*> >& nonIndexedLoopRelations) {
            using namespace std;

            bool compatible = true;

            /**
             * Must order equation pairs property according to the highest 
             * number of shared operations
             */
            map<size_t, map<UniqueEquationPair<Base>, TotalOps2validDepsType*> > totalOp2eq;

            typename set<EquationPattern<Base>*>::const_iterator ite1;
            typename set<EquationPattern<Base>*>::const_iterator ite2;

            for (ite1 = loop1->equations.begin(); ite1 != loop1->equations.end(); ++ite1) {
                EquationPattern<Base>* eq1 = *ite1;

                for (ite2 = loop2->equations.begin(); ite2 != loop2->equations.end(); ++ite2) {
                    EquationPattern<Base>* eq2 = *ite2;

                    UniqueEquationPair<Base> eqRel(eq1, eq2);

                    typename Eq2totalOps2validDepsType::const_iterator eqSharedit = eq2totalOps2validDeps.find(eqRel);
                    if (eqSharedit == eq2totalOps2validDeps.end())
                        continue; // nothing is shared between eq1 and eq2


                    size_t maxOps = eqSharedit->second->rbegin()->first;
                    totalOp2eq[maxOps][eqRel] = eqSharedit->second;
                }
            }

            typename map<size_t, map<UniqueEquationPair<Base>, TotalOps2validDepsType*> >::const_reverse_iterator itr;
            for (itr = totalOp2eq.rbegin(); itr != totalOp2eq.rend(); ++itr) {
                // loop shared operation count

                typename map<UniqueEquationPair<Base>, TotalOps2validDepsType*>::const_iterator itEq;
                for (itEq = itr->second.begin(); itEq != itr->second.end(); ++itEq) {
                    EquationPattern<Base>* eq1 = itEq->first.eq1;
                    EquationPattern<Base>* eq2 = itEq->first.eq2;
                    TotalOps2validDepsType& totalOps2validDeps = *itEq->second;

                    /***************************************************
                     * attempt to combine dependents which share the 
                     * highest number of operations first
                     **************************************************/
                    typename map<size_t, map<DepPairType, const map<OperationNode<Base>*, Indexed2OpCountType>* > >::const_reverse_iterator itOp2Dep2Shared;
                    for (itOp2Dep2Shared = totalOps2validDeps.rbegin(); itOp2Dep2Shared != totalOps2validDeps.rend(); ++itOp2Dep2Shared) {
#ifdef CPPADCG_PRINT_DEBUG
                        std::cout << "    operation count: " << itOp2Dep2Shared->first << "  relations: " << itOp2Dep2Shared->second.size() << std::endl;
#endif
                        typename map<DepPairType, const map<OperationNode<Base>*, Indexed2OpCountType>* >::const_iterator itDep2Shared;
                        for (itDep2Shared = itOp2Dep2Shared->second.begin(); itDep2Shared != itOp2Dep2Shared->second.end(); ++itDep2Shared) {
                            DepPairType depRel = itDep2Shared->first;
                            size_t dep1 = depRel.first;
                            size_t dep2 = depRel.second;

                            const map<OperationNode<Base>*, Indexed2OpCountType>& shared = *itDep2Shared->second;
                            /**
                             * this dep1 <-> dep2 is used as a reference to combine
                             * the two equations in the same loop
                             */
                            compatible = findDepRelations(eq1, dep1, eq2, dep2, shared,
                                                          dep2Relations, dependentBlackListRelations, dependentRelations);
                            if (!compatible) break;
                        }

                        loopRelations.clear();

                        if (compatible) {
                            /**
                             * there has to be at least one iteration with all equation patterns
                             */
                            typename set<EquationPattern<Base>*>::const_iterator ite;
                            std::vector<Loop<Base>*> loops(2);
                            loops[0] = loop1;
                            loops[1] = loop2;
                            bool nonIndexedOnly = true;
                            for (size_t l = 0; l < 2; l++) {
                                Loop<Base>* loop = loops[l];
                                for (ite = loop->equations.begin(); ite != loop->equations.end(); ++ite) { // equation
                                    EquationPattern<Base>* eq = *ite;

                                    set<size_t>::const_iterator itd;
                                    for (itd = eq->dependents.begin(); itd != eq->dependents.end(); ++itd) { // dependent
                                        size_t dep = *itd;
                                        if (dep2Relations[dep] != NULL) {
                                            loopRelations.insert(dep2Relations[dep]);
                                            nonIndexedOnly = false;
                                        }
                                    }
                                }
                            }



                            if (nonIndexedOnly) {
                                nonIndexedLoopRelations.push_back(std::make_pair(eq1, eq2));
                            } else {
                                // there are shared indexed temporary variables
                                compatible = false;
                                set<set<size_t>*>::const_iterator itit;
                                size_t nNonIndexedRel1 = loop1->getLinkedEquationsByNonIndexedCount();
                                size_t nNonIndexedRel2 = loop2->getLinkedEquationsByNonIndexedCount();
                                size_t requiredSize = loop1->equations.size() + loop2->equations.size() - nNonIndexedRel1 - nNonIndexedRel2;

                                for (itit = loopRelations.begin(); itit != loopRelations.end(); ++itit) {
                                    set<size_t>* relations = *itit;
                                    if (relations->size() == requiredSize) {
                                        compatible = true;
                                        break;
                                    }
                                }
#ifdef CPPADCG_PRINT_DEBUG
                                if (compatible) {
                                    std::cout << "    loopRelations:";
                                    print(loopRelations);
                                    std::cout << std::endl;
                                }
#endif
                            }
                        }

                        if (!compatible) break;
                    }

                    if (!compatible) {
                        incompatible_[eq1].insert(eq2);
                        incompatible_[eq2].insert(eq1);
                        break;
                    } else {
                        indexedLoopRelations.insert(eq1);
                        indexedLoopRelations.insert(eq2);
                    }
                }

                if (!compatible) break;
            }

            return compatible;
        }

        /**
         * Determines the relations between the dependents of two equation 
         * patterns using two dependents as reference and the shared variables
         * among the two patterns
         * 
         * @return true if these equation patterns are compatible 
         *         (if they can be in the same loop)
         */
        bool findDepRelations(EquationPattern<Base>* eq1,
                              size_t dep1,
                              EquationPattern<Base>* eq2,
                              size_t dep2,
                              const std::map<OperationNode<Base>*, Indexed2OpCountType>& sharedNodes,
                              std::vector<std::set<size_t>* >& dep2Relations,
                              std::map<size_t, std::set<size_t> >& dependentBlackListRelations,
                              SmartSetPointer<std::set<size_t> >& dependentRelations) {
            using namespace std;

            typename map<OperationNode<Base>*, Indexed2OpCountType>::const_iterator itShared;
            for (itShared = sharedNodes.begin(); itShared != sharedNodes.end(); ++itShared) {
                OperationNode<Base>* sharedNode = itShared->first;

                // checks independents
                bool compatible = canCombineEquations(*eq1, dep1, *eq2, dep2, *sharedNode,
                                                      dep2Relations, dependentBlackListRelations, dependentRelations);

                if (!compatible) return false;
            }

            return true;
        }

        void groupByLoopEqOp(EquationPattern<Base>* eq,
                             std::map<Loop<Base>*, std::map<EquationPattern<Base>*, std::map<size_t, std::pair<OperationNode<Base>*, bool> > > >& loopSharedTemps,
                             const std::map<OperationNode<Base>*, std::set<size_t> >& opShared,
                             bool indexed) {
            using namespace std;

            typename set<OperationNode<Base>*>::const_iterator itShared;
            for (itShared = opShared.begin(); itShared != opShared.end(); ++itShared) {
                OperationNode<Base>* shared = *itShared;
                const set<size_t>& deps = id2Deps[shared->getVariableID()];

                for (set<size_t>::const_iterator itDeps = deps.begin(); itDeps != deps.end(); ++itDeps) {
                    size_t dep = *itDeps;
                    EquationPattern<Base>* otherEq = dep2Equation_.at(dep);
                    if (eq != otherEq) {
                        Loop<Base>* loop = equation2Loop_.at(otherEq);
                        // the original ID (saved in evaluation order) is used to sort shared variables
                        // to ensure reproducibility between different runs
                        loopSharedTemps[loop][otherEq][shared->getEvaluationOrder()] = std::make_pair(shared, indexed);
                    }
                }
            }
        }

        /**
         * Creates a new tape for the model without the equations in the loops
         * and with some extra dependents for the temporary variables used by
         * loops.
         * 
         * @return The new tape without loop equations
         */
        virtual LoopFreeModel<Base>* createNewTape() {
            CodeHandler<Base>& origHandler = *independents_[0].getCodeHandler();

            size_t m = dependents_.size();
            std::vector<bool> inLoop(m, false);
            size_t eqInLoopCount = 0;

            /**
             * Create the new tape
             */
            size_t l_size = loops_.size();

            for (size_t l = 0; l < l_size; l++) {
                Loop<Base>* loop = loops_[l];
                LoopModel<Base>* loopModel = loop->getModel();

                /**
                 * determine which equations belong to loops
                 */
                const std::vector<std::vector<LoopPosition> >& ldeps = loopModel->getDependentIndexes();
                for (size_t eq = 0; eq < ldeps.size(); eq++) {
                    for (size_t it = 0; it < ldeps[eq].size(); it++) {
                        const LoopPosition& pos = ldeps[eq][it];
                        if (pos.original != std::numeric_limits<size_t>::max()) {// some equations are not present in all iteration
                            inLoop[pos.original] = true;
                            eqInLoopCount++;
                        }
                    }
                }
            }

            /**
             * create a new smaller tape 
             */
            size_t nonLoopEq = m - eqInLoopCount;
            std::vector<CGBase> nonLoopDeps(nonLoopEq + origTemp2Index_.size());

            if (nonLoopDeps.size() == 0)
                return NULL; // there are no equations outside the loops

            /**
             * Place the dependents that do not belong to a loop
             */
            size_t inl = 0;
            std::vector<size_t> depTape2Orig(nonLoopEq);
            if (eqInLoopCount < m) {
                for (size_t i = 0; i < inLoop.size(); i++) {
                    if (!inLoop[i]) {
                        depTape2Orig[inl] = i;
                        nonLoopDeps[inl++] = dependents_[i];
                    }
                }
            }
            CPPADCG_ASSERT_UNKNOWN(inl == nonLoopEq);

            /**
             * Place new dependents for the temporary variables used by the loops
             */
            typename std::map<OperationNode<Base>*, size_t>::const_iterator itTmp;
            for (itTmp = origTemp2Index_.begin(); itTmp != origTemp2Index_.end(); ++itTmp) {
                size_t k = itTmp->second;
                nonLoopDeps[nonLoopEq + k] = origHandler.createCG(Argument<Base>(*itTmp->first));
            }

            /**
             * Generate the new tape by going again through the operations 
             */
            Evaluator<Base, CGBase> evaluator(origHandler, nonLoopDeps);

            // set atomic functions
            const std::map<size_t, CGAbstractAtomicFun<Base>* >& atomicsOrig = origHandler.getAtomicFunctions();
            std::map<size_t, atomic_base<CGBase>* > atomics;
            atomics.insert(atomicsOrig.begin(), atomicsOrig.end());
            evaluator.addAtomicFunctions(atomics);

            std::vector<AD<CGBase> > x(independents_.size());
            for (size_t j = 0; j < x.size(); j++) {
                if (independents_[j].isValueDefined())
                    x[j] = independents_[j].getValue();
            }

            CppAD::Independent(x);
            std::vector<AD<CGBase> > y = evaluator.evaluate(x);

            std::unique_ptr<ADFun<CGBase> > tapeNoLoops(new ADFun<CGBase>());
            tapeNoLoops->Dependent(y);

            return new LoopFreeModel<Base>(tapeNoLoops.release(), depTape2Orig);
        }

        std::vector<EquationPattern<Base>*> findRelatedVariables() throw (CGException) {
            eqCurr_ = NULL;
            color_ = 1; // used to mark visited nodes

            size_t rSize = relatedDepCandidates_.size();
            for (size_t r = 0; r < rSize; r++) {
                const std::set<size_t>& candidates = relatedDepCandidates_[r];
                std::set<size_t> used;

                eqCurr_ = NULL;

                std::set<size_t>::const_iterator itRef;
                for (itRef = candidates.begin(); itRef != candidates.end(); ++itRef) {
                    size_t iDepRef = *itRef;

                    // check if it has already been used
                    if (used.find(iDepRef) != used.end()) {
                        continue;
                    }

                    if (eqCurr_ == NULL || used.size() > 0) {
                        eqCurr_ = new EquationPattern<Base>(dependents_[iDepRef], iDepRef);
                        equations_.push_back(eqCurr_);
                    }

                    std::set<size_t>::const_iterator it = itRef;
                    for (++it; it != candidates.end(); ++it) {
                        size_t iDep = *it;
                        // check if it has already been used
                        if (used.find(iDep) != used.end()) {
                            continue;
                        }

                        if (eqCurr_->testAdd(iDep, dependents_[iDep], color_)) {
                            used.insert(iDep);
                        }
                    }

                    if (eqCurr_->dependents.size() == 1) {
                        // nothing found :(
                        delete eqCurr_;
                        eqCurr_ = NULL;
                        equations_.pop_back();
                    }
                }
            }

            /**
             * Determine the independents that don't change from iteration to
             * iteration
             */
            for (size_t eq = 0; eq < equations_.size(); eq++) {
                equations_[eq]->detectNonIndexedIndependents();
            }

            return equations_;
        }

        /**
         * Finds nodes which can be shared with other equation patterns
         * 
         * @param value The CG object to visit
         * @param depIndex The index of the dependent variable
         * @return true if this operation is indexed
         */
        inline bool findSharedTemporaries(const CG<Base>& value,
                                          size_t depIndex) {
            OperationNode<Base>* depNode = value.getOperationNode();
            size_t opCount = 0;
            if (findSharedTemporaries(depNode, depIndex, opCount)) {
                depNode->setColor(color_);
                return true;
            }
            return false;
        }

        /**
         * Finds nodes which can be shared with other equation patterns
         * 
         * @param node The node to visit
         * @param depIndex The index of the dependent variable
         * @param opCount The number of operations which must be performed to 
         *                generate this node
         * @return true if this operation is indexed (for the current equation pattern)
         */
        inline bool findSharedTemporaries(OperationNode<Base>* node,
                                          size_t depIndex,
                                          size_t& opCount) {
            if (node == NULL)
                return false; // nothing to do

            if (node->getUsageCount() > 0) {
                opCount++; // this operation
                return node->getColor() == color_;
            }

            node->increaseUsageCount();

            bool indexedOperation = false;

            size_t localOpCount = 1;
            const std::vector<Argument<Base> >& args = node->getArguments();
            size_t arg_size = args.size();
            for (size_t a = 0; a < arg_size; a++) {
                OperationNode<Base>*argOp = args[a].getOperation();
                if (argOp != NULL) {
                    if (argOp->getOperationType() != CGInvOp) {
                        indexedOperation |= findSharedTemporaries(argOp, depIndex, localOpCount);
                    } else {
                        indexedOperation |= !eqCurr_->containsConstantIndependent(node, a);
                    }
                }
            }

            opCount += localOpCount;

            if (indexedOperation) {
                node->setColor(color_); // mark this operation as being indexed
            } else {
                node->setColor(0); // mark this operation as being not-indexed
            }

            size_t id = node->getVariableID();
            std::set<size_t>& deps = id2Deps[id];

            if (deps.size() > 1 && node->getOperationType() != CGInvOp) {
                /**
                 * Temporary variable
                 */
                std::set<size_t>::const_iterator it;
                for (it = deps.begin(); it != deps.end(); ++it) {
                    size_t otherDep = *it;

                    EquationPattern<Base>* otherEquation = dep2Equation_.at(otherDep);
                    if (otherEquation != eqCurr_) {
                        /**
                         * temporary variable shared with a different loop
                         */
                        UniqueEquationPair<Base> eqPair(eqCurr_, otherEquation);
                        Dep1Dep2SharedType& relation = equationShared_[eqPair];

                        std::map<OperationNode<Base>*, Indexed2OpCountType>* reldepdep;
                        if (eqPair.eq1 == eqCurr_)
                            reldepdep = &relation[depIndex][otherDep];
                        else
                            reldepdep = &relation[otherDep][depIndex];

                        INDEXED_OPERATION_TYPE expected = indexedOperation ? INDEXED_OPERATION_TYPE_INDEXED : INDEXED_OPERATION_TYPE_NONINDEXED;
                        typename std::map<OperationNode<Base>*, Indexed2OpCountType>::iterator itIndexedType = reldepdep->find(node);
                        if (itIndexedType == reldepdep->end()) {
                            (*reldepdep)[node] = Indexed2OpCountType(expected, localOpCount);
                        } else if (itIndexedType->second.first != expected) {
                            itIndexedType->second.first = INDEXED_OPERATION_TYPE_BOTH;
                        }

                        break;
                    }
                }
            }

            return indexedOperation;
        }

        /**
         * Marks a node (and all other nodes used by it) as being used by a
         * given dependent variable.
         * 
         * @param node The node being visited
         * @param dep Dependent variable index
         */
        inline void markOperationsWithDependent(const OperationNode<Base>* node,
                                                size_t dep) {
            if (node == NULL || node->getOperationType() == CGInvOp)
                return; // nothing to do

            size_t id = node->getVariableID();

            std::set<size_t>& deps = id2Deps[id];

            if (deps.size() == 0) {
                deps.insert(dep); // here for the first time 
            } else {
                std::pair < std::set<size_t>::iterator, bool> added = deps.insert(dep);
                if (!added.second) {
                    return; // already been here
                }
            }

            const std::vector<Argument<Base> >& args = node->getArguments();
            size_t arg_size = args.size();
            for (size_t i = 0; i < arg_size; i++) {
                markOperationsWithDependent(args[i].getOperation(), dep);
            }
        }

        void assignIds() {
            idCounter_ = 1;

            size_t rSize = relatedDepCandidates_.size();
            for (size_t r = 0; r < rSize; r++) {
                const std::set<size_t>& candidates = relatedDepCandidates_[r];

                std::set<size_t>::const_iterator it;
                for (it = candidates.begin(); it != candidates.end(); ++it) {
                    assignIds(dependents_[*it].getOperationNode());
                }
            }
        }

        void assignIds(OperationNode<Base>* node) {
            if (node == NULL || node->getVariableID() > 0)
                return;

            node->setVariableID(idCounter_);
            node->setEvaluationOrder(idCounter_);
            idCounter_++;

            const std::vector<Argument<Base> >& args = node->getArguments();
            size_t arg_size = args.size();
            for (size_t i = 0; i < arg_size; i++) {
                assignIds(args[i].getOperation());
            }
        }

        void resetHandlerCounters() {
            size_t rSize = relatedDepCandidates_.size();
            for (size_t r = 0; r < rSize; r++) {
                const std::set<size_t>& candidates = relatedDepCandidates_[r];

                std::set<size_t>::const_iterator it;
                for (it = candidates.begin(); it != candidates.end(); ++it) {
                    resetHandlerCounters(dependents_[*it].getOperationNode());
                }
            }
        }

        static void resetHandlerCounters(OperationNode<Base>* node) {
            if (node == NULL || node->getVariableID() == 0 || node->getEvaluationOrder() == 0)
                return;

            node->resetHandlerCounters();

            const std::vector<Argument<Base> >& args = node->getArguments();
            size_t arg_size = args.size();
            for (size_t i = 0; i < arg_size; i++) {
                resetHandlerCounters(args[i].getOperation());
            }
        }

        static bool find(Loop<Base>* loop1, Loop<Base>* loop2,
                         const std::map<EquationPattern<Base>*, std::set<EquationPattern<Base>*> >& blackList) {
            typename std::set<EquationPattern<Base>*>::const_iterator iteq1;
            for (iteq1 = loop1->equations.begin(); iteq1 != loop1->equations.end(); ++iteq1) {

                typename std::map<EquationPattern<Base>*, std::set<EquationPattern<Base>* > >::const_iterator itBlack;
                itBlack = blackList.find(*iteq1);
                if (itBlack != blackList.end()) {

                    typename std::set<EquationPattern<Base>*>::const_iterator iteq2;
                    for (iteq2 = loop2->equations.begin(); iteq2 != loop2->equations.end(); ++iteq2) {
                        if (itBlack->second.find(*iteq2) != itBlack->second.end()) {
                            return true; // found
                        }
                    }
                }
            }

            return false;
        }

        template<class T>
        static inline bool contains(const std::map<T, std::set<T> >& map, T eq1, T eq2) {
            typename std::map<T, std::set<T> >::const_iterator itb1;
            itb1 = map.find(eq1);
            if (itb1 != map.end()) {
                if (itb1->second.find(eq2) != itb1->second.end()) {
                    return true;
                }
            }
            return false;
        }

        bool canCombineEquations(const EquationPattern<Base>& eq1,
                                 size_t dep1,
                                 const EquationPattern<Base>& eq2,
                                 size_t dep2,
                                 OperationNode<Base>& sharedTemp,
                                 std::vector<std::set<size_t>* >& dep2Relations,
                                 std::map<size_t, std::set<size_t> >& dependentBlackListRelations,
                                 SmartSetPointer<std::set<size_t> >& dependentRelations) {
            using namespace std;

            // must have indexed independents at the same locations in all equations
            const set<const OperationNode<Base>*> opWithIndepArgs = EquationPattern<Base>::findOperationsUsingIndependents(sharedTemp);
            EquationPattern<Base>::clearUsageCount(&sharedTemp); // must reset usage count

            // must have indexed independents at the same locations in both equations
            typename set<const OperationNode<Base>*>::const_iterator itOp;
            for (itOp = opWithIndepArgs.begin(); itOp != opWithIndepArgs.end(); ++itOp) {
                const OperationNode<Base>* op = *itOp;

                // get indexed independent variable information
                // - equation 1
                typename map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::const_iterator indexed1It;
                OperationNode<Base>* op1 = eq1.operationEO2Reference.at(dep1).at(op); // convert to the reference of equation 1
                indexed1It = eq1.indexedOpIndep.op2Arguments.find(op1);

                // - equation 2
                typename map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::const_iterator indexed2It;
                OperationNode<Base>* op2 = eq2.operationEO2Reference.at(dep2).at(op); // convert to the reference of equation 2
                indexed2It = eq2.indexedOpIndep.op2Arguments.find(op2);

                /**
                 * Compare the iterations where this operation is used in both equation patterns
                 */
                if (indexed1It == eq1.indexedOpIndep.op2Arguments.end()) {
                    if (indexed2It != eq2.indexedOpIndep.op2Arguments.end()) {
                        return false; // indexed in one equation but non-indexed in the other
                    }
                } else {
                    if (indexed2It == eq2.indexedOpIndep.op2Arguments.end()) {
                        return false; // indexed in one equation but non-indexed in the other
                    }

                    /**
                     * Indexed path in both equations
                     */
                    const OperationIndexedIndependents<Base>& indexed1Ops = indexed1It->second;
                    const OperationIndexedIndependents<Base>& indexed2Ops = indexed2It->second;

                    size_t a1Size = indexed1Ops.arg2Independents.size();
                    if (a1Size != indexed2Ops.arg2Independents.size()) { // there must be the same number of arguments
                        return false;
                    }

                    for (size_t a = 0; a < a1Size; a++) {
                        const map<size_t, const OperationNode<Base>*>& eq1Dep2Indep = indexed1Ops.arg2Independents[a];
                        const map<size_t, const OperationNode<Base>*>& eq2Dep2Indep = indexed2Ops.arg2Independents[a];

                        if (eq1Dep2Indep.empty() != eq2Dep2Indep.empty())
                            return false; // one is indexed and the other is non-indexed

                        // it has to be possible to match dependents from the two equation patterns

                        if (eq1Dep2Indep.empty()) {
                            continue; // not indexed
                        }

                        // indexed independent variable

                        // invert eq1Dep2Indep into eq1Indep2Dep
                        typedef map<const OperationNode<Base>*, size_t, IndependentNodeSorter<Base> > MapIndep2Dep;
                        MapIndep2Dep eq1Indep2Dep;
                        typename MapIndep2Dep::iterator hint = eq1Indep2Dep.begin();
                        typename map<size_t, const OperationNode<Base>*>::const_iterator d2i;
                        for (d2i = eq1Dep2Indep.begin(); d2i != eq1Dep2Indep.end(); ++d2i) {
                            hint = eq1Indep2Dep.insert(hint, std::make_pair(d2i->second, d2i->first));
                            hint++; // assume that the relation dep<->indep is always ascending
                        }

                        typename map<const OperationNode<Base>*, size_t>::const_iterator itHint = eq1Indep2Dep.begin();

                        // check all iterations/dependents
                        for (d2i = eq2Dep2Indep.begin(); d2i != eq2Dep2Indep.end(); ++d2i) {
                            size_t dep2 = d2i->first;
                            const OperationNode<Base>* indep = d2i->second;
                            typename map<const OperationNode<Base>*, size_t>::const_iterator it;
                            if (itHint->first == indep) {
                                /**
                                 * assumes that the both relations 
                                 * (indep<->dep1 and dep2<->indep) are in
                                 * ascending order which is commonly true
                                 */
                                it = itHint;
                                itHint++;
                            } else {
                                it = eq1Indep2Dep.find(indep);
                            }

                            if (it != eq1Indep2Dep.end()) {
                                size_t dep1 = it->second;

                                // check if this relation was previous excluded
                                std::map<size_t, set<size_t> >::const_iterator itBlackL = dependentBlackListRelations.find(dep1);
                                if (itBlackL != dependentBlackListRelations.end() && itBlackL->second.find(dep2) != itBlackL->second.end()) {
                                    return false; // these dependents cannot be in the same iteration
                                }

                                bool related = makeDependentRelation(eq1, dep1, eq2, dep2,
                                                                     dep2Relations, dependentRelations);
                                if (!related)
                                    return false;

                            } else {
                                // equation pattern 1 does not have any iteration with indep from dep2

                                // there is no need to have the same number of iterations in both equations!
                                // but remember that these dependents cannot be in the same iteration from now on
                                dependentBlackListRelations[dep2].insert(eq1.dependents.begin(), eq1.dependents.end());
                            }
                        }

                    }

                }

            }

            return true;
        }

        bool isNonIndexed(const EquationPattern<Base>& eq2,
                          size_t dep2,
                          OperationNode<Base>& sharedTemp) {
            using namespace std;


            // must have indexed independents at the same locations in all equations
            const set<const OperationNode<Base>*> opWithIndepArgs = EquationPattern<Base>::findOperationsUsingIndependents(sharedTemp);

            typename set<const OperationNode<Base>*>::const_iterator itOp;
            for (itOp = opWithIndepArgs.begin(); itOp != opWithIndepArgs.end(); ++itOp) {
                const OperationNode<Base>* op = *itOp;

                // get indexed independent variable information
                // - equation 2
                typename map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::const_iterator indexed2It;
                OperationNode<Base>* op2 = eq2.operationEO2Reference.at(dep2).at(op); // convert to the reference of equation 2
                indexed2It = eq2.indexedOpIndep.op2Arguments.find(op2);

                if (indexed2It != eq2.indexedOpIndep.op2Arguments.end()) {
                    return false; // indexed in one equation but non-indexed in the other
                }
            }

            return true;
        }

        bool makeDependentRelation(const EquationPattern<Base>& eq1,
                                   size_t dep1,
                                   const EquationPattern<Base>& eq2,
                                   size_t dep2,
                                   std::vector<std::set<size_t>* >& dep2Relations,
                                   SmartSetPointer<std::set<size_t> >& dependentRelations) {
            using namespace std;

            std::set<size_t>::const_iterator it;
            set<size_t>* related1 = dep2Relations[dep1];
            set<size_t>* related2 = dep2Relations[dep2];

            // check if relations were established with a different dependent from the same equation pattern
            if (related1 != NULL) {
                // dependent 1 already in a relation set

                if (related2 != NULL) {
                    // both dependents belong to previously existing relations sets

                    if (related1 == related2)
                        return true; // already done

                    // relations must be merged (if possible)!
                    // merge related2 into related1
                    bool canMerge = true;
                    set<size_t>::const_iterator itr2;
                    for (itr2 = related2->begin(); itr2 != related2->end(); ++itr2) {
                        size_t dep3 = *itr2;
                        const EquationPattern<Base>& eq3 = *dep2Equation_.at(dep3);
                        // make sure no other dependent from the same equation pattern was already in this relation set
                        for (it = eq3.dependents.begin(); it != eq3.dependents.end(); ++it) {
                            if (*it != dep3 && related1->find(*it) != related1->end()) {
                                canMerge = false; // relation with a dependent from a different iteration!
                                break;
                                //return false; 
                            }
                        }

                        if (!canMerge)
                            break;
                    }

                    if (canMerge) {
                        for (itr2 = related2->begin(); itr2 != related2->end(); ++itr2) {
                            size_t dep3 = *itr2;
                            related1->insert(dep3);
                            dep2Relations[dep3] = related1;
                        }

                        dependentRelations.erase(related2);
                        delete related2;
                    }
                    /**
                     * when it is not possible to merge due to a dependent
                     * variable from another iteration already belonging to a 
                     * dependent variable relation (iteration) the new relation
                     * is not added and as a consequence there will be some 
                     * repeated operations
                     */

                } else {
                    if (related1->find(dep2) == related1->end()) {
                        // make sure no other dependent from the same equation pattern was already in this relation set
                        bool canMerge = true;
                        for (it = eq2.dependents.begin(); it != eq2.dependents.end(); ++it) {
                            if (*it != dep2 && related1->find(*it) != related1->end()) {
                                canMerge = false; // relation with a dependent from a different iteration!
                                break;
                            }
                        }

                        if (canMerge) {
                            related1->insert(dep2);
                            dep2Relations[dep2] = related1;
                        }
                        /**
                         * when it is not possible to merge due to a dependent
                         * variable from another iteration already belonging to a 
                         * dependent variable relation (iteration) the new relation
                         * is not added and as a consequence there will be some 
                         * repeated operations
                         */
                    }
                }

            } else if (related2 != NULL) {
                // dependent 2 already in a relation set

                // make sure no other dependent from the same equation pattern was already in this relation set
                bool canMerge = true;
                for (it = eq1.dependents.begin(); it != eq1.dependents.end(); ++it) {
                    if (*it != dep1 && related2->find(*it) != related2->end()) {
                        canMerge = false; // relation with a dependent from a different iteration!
                        break;
                        //return false;
                    }
                }

                if (canMerge) {
                    related2->insert(dep1);
                    dep2Relations[dep1] = related2;
                }
                /**
                 * when it is not possible to merge due to a dependent
                 * variable from another iteration already belonging to a 
                 * dependent variable relation (iteration) the new relation
                 * is not added and as a consequence there will be some 
                 * repeated operations
                 */


            } else {
                // dependent 1 and dependent 2 not in any relation set
                set<size_t>* related = new std::set<size_t>();
                dependentRelations.insert(related);
                related->insert(dep1);
                related->insert(dep2);
                dep2Relations[dep1] = related;
                dep2Relations[dep2] = related;
            }

            return true;
        }

    };

}

#endif

