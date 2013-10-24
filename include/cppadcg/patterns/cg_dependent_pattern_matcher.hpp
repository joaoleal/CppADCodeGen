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

namespace CppAD {

    /**
     * Finds common patterns in operation graphs
     */
    template<class Base>
    class DependentPatternMatcher {
    public:
        typedef CG<Base> CGBase;
    private:
        const std::vector<std::set<size_t> >& relatedDepCandidates_;
        const std::vector<CGBase>& dependents_;
        std::vector<CGBase>& independents_;
        std::vector<EquationPattern<Base>*> equations_;
        EquationPattern<Base>* eqCurr_;
        std::map<size_t, EquationPattern<Base>*> dep2Equation_;
        std::map<EquationPattern<Base>*, Loop<Base>*> equation2Loop_;
        std::vector<Loop<Base>*> loops_;
        std::set<const OperationNode<Base>*> sharedTemps_; // shared temporaries
        /**
         * maps the orignal model nodes used as temporary non-indexed variables
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
            assert(independents_.size() > 0);
            assert(independents_[0].getCodeHandler() != NULL);
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

                std::set<size_t>::const_iterator depIt;
                for (depIt = eq->dependents.begin(); depIt != eq->dependents.end(); ++depIt) {
                    dep2Equation_[*depIt] = eq;
                }
            }

            loops_.reserve(eq_size);

            std::map<EquationPattern<Base>*, std::set<Loop<Base>* > > blackList;
            std::map<EquationPattern<Base>*, std::set<EquationPattern<Base>* > > blackListEq;

            /**
             * Combine related equations in the same loops
             * (equations that share temporary variables)
             */
            color_++; // this is the color used to mark indexed nodes (must be higher than any previously used color)

            for (size_t e = 0; e < eq_size; e++) {
                EquationPattern<Base>* eq = equations_[e];
                eqCurr_ = eq;
                std::set<size_t>::const_iterator depIt;

                for (depIt = eq->dependents.begin(); depIt != eq->dependents.end(); ++depIt) {
                    OperationNode<Base>* node = dependents_[*depIt].getOperationNode();
                    // will define the dependents associated with each operation
                    markOperationsWithDependent(node, *depIt);
                }

                /**
                 * Find shared values in the indexed path
                 */
                sharedTemps_.clear();
                for (depIt = eq->dependents.begin(); depIt != eq->dependents.end(); ++depIt) {
                    findSharedTemporaries(dependents_[*depIt]);
                }

                /**
                 * clean-up
                 */
                for (depIt = eq->dependents.begin(); depIt != eq->dependents.end(); ++depIt) {
                    OperationNode<Base>* node = dependents_[*depIt].getOperationNode();
                    // must uncolor
                    EquationPattern<Base>::uncolor(node);
                    // must reset usage count
                    EquationPattern<Base>::clearUsageCount(node);
                }

                if (sharedTemps_.empty()) {
                    /**
                     * new loop (with no relations with other equations)
                     */
                    Loop<Base>* loop = new Loop<Base>(*eq);
                    loops_.push_back(loop);
                    equation2Loop_[eq] = loop;
                } else {

                    Loop<Base>* currLoop = NULL;

                    std::map<Loop<Base>*, std::vector<std::set<size_t> > > newLoopRelations;
                    std::map < Loop<Base>*, bool> canAdd2Loop;
                    if (loops_.size() > 0) {
                        typename std::set<const OperationNode<Base>*>::const_iterator itShared;
                        for (itShared = sharedTemps_.begin(); itShared != sharedTemps_.end(); ++itShared) {
                            const OperationNode<Base>* shared = *itShared;
                            const std::set<size_t>& deps = id2Deps[shared->getVariableID()];
                            std::map<Loop<Base>*, std::set<size_t> > loopRelations4Shared; // relationships due to this shared variable

                            /**
                             * count the number of times the shared variable was 
                             * used by each equation pattern
                             */
                            std::map<EquationPattern<Base>*, std::set<size_t> > eq2depWithShared;
                            for (std::set<size_t>::const_iterator itDeps = deps.begin(); itDeps != deps.end(); ++itDeps) {
                                size_t dep = *itDeps;
                                EquationPattern<Base>* otherEq = dep2Equation_.at(dep);
                                eq2depWithShared[otherEq].insert(dep);
                            }

                            size_t timesShared = eq2depWithShared[eq].size();

                            typename std::map<EquationPattern<Base>*, std::set<size_t> >::const_iterator it;
                            for (it = eq2depWithShared.begin(); it != eq2depWithShared.end(); ++it) {
                                EquationPattern<Base>* otherEq = it->first;
                                const std::set<size_t>& deps = it->second;

                                if (otherEq == eq)
                                    continue;

                                bool compatibleEqs = otherEq->dependents.size() == eq->dependents.size() && // must have the appropriate number of iterations
                                        deps.size() == timesShared; // must show up the same number of times otherwise it is an indexed variable in one equation and non-indexed in another

                                if (!compatibleEqs) {
                                    blackListEq[eq].insert(otherEq);
                                    blackListEq[otherEq].insert(eq);
                                }

                                typename std::map<EquationPattern<Base>*, Loop<Base>*>::const_iterator ite2l = equation2Loop_.find(otherEq);
                                if (ite2l != equation2Loop_.end()) {
                                    // equation already belongs to a loop
                                    Loop<Base>* loop = ite2l->second;

                                    canAdd2Loop[loop] = compatibleEqs;

                                    if (compatibleEqs && timesShared == 1) {
                                        // must have only one dependent per equation otherwise it is not an indexed shared variable
                                        loopRelations4Shared[loop].insert(*deps.begin());
                                    }
                                }
                            }

                            size_t expectedCount = eq2depWithShared.at(eq).size();

                            for (it = eq2depWithShared.begin(); it != eq2depWithShared.end(); ++it) {
                                EquationPattern<Base>* otherEq = it->first;
                                typename std::map<EquationPattern<Base>*, Loop<Base>*>::const_iterator ite2l = equation2Loop_.find(otherEq);
                                if (ite2l != equation2Loop_.end()) {
                                    Loop<Base>* loop = ite2l->second;
                                    if (canAdd2Loop[loop]) {
                                        if (it->second.size() != expectedCount) {
                                            // different number of times and therefore it is only shared by some indexes/dependents
                                            canAdd2Loop[loop] = false;
                                        } else {
                                            /////////// findOperationsUsingIndependents() called too many times for the same shared variables//////////////////////////////
                                            // check if it is indexed in same locations 
                                            canAdd2Loop[loop] |= Loop<Base>::canCombineEquations(*eq, *otherEq, *shared); // checks independents
                                        }
                                    }
                                }
                            }

                            typename std::map<Loop<Base>*, std::set<size_t> >::const_iterator itlr;
                            for (itlr = loopRelations4Shared.begin(); itlr != loopRelations4Shared.end(); ++itlr) {
                                Loop<Base>* loop = itlr->first;
                                if (canAdd2Loop[loop]) {
                                    Loop<Base>::combineDependentRelations(itlr->second, newLoopRelations[loop]);
                                }
                            }
                        }

                        typename std::map < Loop<Base>*, bool>::const_iterator itLoop;
                        for (itLoop = canAdd2Loop.begin(); itLoop != canAdd2Loop.end(); ++itLoop) {
                            Loop<Base>* loop = itLoop->first;

                            if (!itLoop->second) {
                                blackList[eq].insert(loop);
                                continue; // cannot add the equation to this loop
                            }

                            if (currLoop == NULL) {
                                // add to an existing loop
                                const std::vector<std::set<size_t> >& newDepRel = newLoopRelations[loop];
                                bool ok = loop->isCompatibleWithDepRelationships(newDepRel, dep2Equation_); // checks dependents
                                if (ok) {
                                    currLoop = loop;
                                    currLoop->addEquation(*eq, newDepRel);
                                    equation2Loop_[eq] = loop;
                                } else {
                                    blackList[eq].insert(loop);
                                }
                            } else {
                                // had already been added to a loop
                                // try to merge the two loops if possible
                                bool ok = loop->isCompatibleWithDepRelationships(currLoop->linkedDependents, dep2Equation_); // checks dependents
                                if (ok) {
                                    // update the loop of the equations
                                    typename std::set<EquationPattern<Base>*>::const_iterator itle;
                                    for (itle = loop->equations.begin(); itle != loop->equations.end(); ++itle) {
                                        equation2Loop_[*itle] = currLoop;
                                    }
                                    currLoop->merge(*loop);

                                    typename std::vector<Loop<Base>*>::iterator it = std::find(loops_.begin(), loops_.end(), loop);
                                    assert(it != loops_.end());
                                    loops_.erase(it);
                                    delete loop;
                                } else {
                                    blackList[eq].insert(loop);
                                }
                            }
                        }
                    }

                    if (currLoop == NULL) {
                        // could not add the equation to any existing loop
                        Loop<Base>* loop = new Loop<Base>(*eq);
                        loops_.push_back(loop);
                        equation2Loop_[eq] = loop;
                    }

                }

            }

            /**
             * Attempt to combine unrelated loops
             */
            if (!loops_.empty()) {
                for (size_t l1 = 0; l1 < loops_.size() - 1; l1++) {
                    Loop<Base>* loop1 = loops_[l1];
                    for (size_t l2 = l1 + 1; l2 < loops_.size();) {
                        Loop<Base>* loop2 = loops_[l2];

                        bool canMerge = loop1->iterationCount == loop2->iterationCount;
                        if (canMerge) {
                            // check if there are equations in the blacklist
                            canMerge = !find(loop1, loop2, blackList) && !find(loop2, loop1, blackList);
                            if (canMerge) {
                                typename std::set<EquationPattern<Base>*>::const_iterator itl1, itl2;
                                for (itl1 = loop1->equations.begin(); itl1 != loop1->equations.begin(); ++itl1) {
                                    EquationPattern<Base>* eq1 = *itl1;
                                    for (itl2 = loop2->equations.begin(); itl2 != loop2->equations.begin(); ++itl2) {
                                        EquationPattern<Base>* eq2 = *itl2;
                                        if (contains(blackListEq, eq1, eq2)) {
                                            canMerge = false;
                                            break;
                                        }
                                    }
                                    if (!canMerge)
                                        break;
                                }
                            }
                        }

                        if (canMerge) {
                            loop1->merge(*loop2);
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
                        inLoop[pos.original] = true;
                    }
                }
                eqInLoopCount += ldeps.size() * loopModel->getIterationCount();
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
            assert(inl == nonLoopEq);

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

            std::auto_ptr<ADFun<CGBase> > tapeNoLoops(new ADFun<CGBase>());
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
         * @param value the CG object to visit
         * @return true if this operation is indexed
         */
        inline bool findSharedTemporaries(const CG<Base>& value) {
            OperationNode<Base>* depNode = value.getOperationNode();
            if (findSharedTemporaries(depNode)) {
                depNode->setColor(color_);
                return true;
            }
            return false;
        }

        /**
         * Finds nodes which can be shared with other equation patterns
         * 
         * @param node the node to visit
         * @return true if this operation is indexed
         */
        inline bool findSharedTemporaries(OperationNode<Base>* node) {
            if (node == NULL)
                return false; // nothing to do

            if (node->getUsageCount() > 0)
                return node->getColor() == color_;

            node->increaseUsageCount();

            bool indexedOperation = false;
            const std::vector<Argument<Base> >& args = node->getArguments();
            size_t arg_size = args.size();
            for (size_t a = 0; a < arg_size; a++) {
                OperationNode<Base>*argOp = args[a].getOperation();
                if (argOp != NULL) {
                    if (argOp->getOperationType() != CGInvOp) {
                        indexedOperation |= findSharedTemporaries(argOp);
                    } else {
                        indexedOperation |= !eqCurr_->containsConstantIndependent(node, a);
                    }
                }
            }

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
                        sharedTemps_.insert(node);
                    }
                }
            }

            return indexedOperation;
        }

        inline void markOperationsWithDependent(OperationNode<Base>* node, size_t dep) {
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
                         std::map<EquationPattern<Base>*, std::set<Loop<Base>* > > blackList) {
            typename std::set<EquationPattern<Base>*>::const_iterator iteq;
            for (iteq = loop1->equations.begin(); iteq != loop1->equations.end(); ++iteq) {
                typename std::map<EquationPattern<Base>*, std::set<Loop<Base>* > >::const_iterator itBlack;
                itBlack = blackList.find(*iteq);
                if (itBlack != blackList.end() &&
                        itBlack->second.find(loop2) != itBlack->second.end()) {
                    return true; // found
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

    };

}

#endif

