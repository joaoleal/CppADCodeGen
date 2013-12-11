#ifndef CPPAD_CG_EQUATION_PATTERN_INCLUDED
#define CPPAD_CG_EQUATION_PATTERN_INCLUDED
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
     * Holds information on which independents are used by which 
     * dependent in an equation pattern
     */
    template<class Base>
    class OperationIndexedIndependents {
    public:
        typedef std::map<size_t, const OperationNode<Base>*> MapDep2Indep_type;
        /**
         * maps the argument index to the several independents used in different
         * equations with the same pattern
         * argument index -> (iteration or dependent using it) -> independent
         */
        std::vector<MapDep2Indep_type> arg2Independents;
    };

    template<class Base>
    class IndexedIndependent {
    public:
        std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base> > op2Arguments;
    public:

        const OperationIndexedIndependents<Base>& arguments(const OperationNode<Base>* operation) const {
            return op2Arguments.at(operation);
        }

        bool isIndexedOperationArgument(const OperationNode<Base>* node, size_t argIndex) const {
            typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::const_iterator itIndexes;
            itIndexes = op2Arguments.find(node);
            if (itIndexes == op2Arguments.end()) {
                return false;
            }
            const OperationIndexedIndependents<Base>& indexedArgs = itIndexes->second;
            return indexedArgs.arg2Independents.size() > argIndex && !indexedArgs.arg2Independents[argIndex].empty();
        }

    };

    /**
     * Group of variables with the same evaluation pattern 
     * (same equation different variables)
     */
    template<class Base>
    class EquationPattern {
    public:
        const CG<Base>& depRef; // dependent reference
        const size_t depRefIndex;
        std::set<size_t> dependents;
        /**
         * maps node ID used by all dependents to the operations of the 
         * reference dependent
         * [dependent index][op] = reference operation
         */
        std::map<size_t, std::map<const OperationNode<Base>*, OperationNode<Base>*> > operationEO2Reference;
        // std::map<size_t, std::vector<OperationNode<Base>*> > operationEO2Reference;
        /**
         * Maps the operations that used an indexed independents as direct
         * arguments
         * (reference operation -> argument indexes -> dependent <-> independents)
         */
        IndexedIndependent<Base> indexedOpIndep;
        /**
         * reference operation -> non indexed argument positions
         */
        std::map<const OperationNode<Base>*, std::set<size_t> > constOperationIndependents;

    private:
        size_t currDep_;
        size_t minColor_;
        size_t cmpColor_;
        /**
         * a unique index
         */
        size_t index_;
    public:

        explicit EquationPattern(const CG<Base>& ref,
                                 size_t iDepRef) :
            depRef(ref),
            depRefIndex(iDepRef) {
            dependents.insert(iDepRef);
        }

        bool testAdd(size_t iDep2, const CG<Base>& dep2, size_t& minColor) throw (CGException) {
            IndexedIndependent<Base> independentsBackup = indexedOpIndep;
            std::map<const OperationNode<Base>*, std::set<size_t> > constOperationIndependentsBackup = constOperationIndependents;
            std::map<size_t, std::map<const OperationNode<Base>*, OperationNode<Base>*> > operation2ReferenceBackup = operationEO2Reference;

            currDep_ = iDep2;
            minColor_ = minColor;
            cmpColor_ = minColor_;

            bool equals = comparePath(depRef, dep2, iDep2);

            minColor = cmpColor_;

            if (equals) {
                dependents.insert(iDep2);

                return true; // matches the reference pattern
            } else {
                // restore
                indexedOpIndep.op2Arguments.swap(independentsBackup.op2Arguments);
                constOperationIndependents.swap(constOperationIndependentsBackup);
                operationEO2Reference.swap(operation2ReferenceBackup);

                return false; // cannot be added
            }
        }

        inline void colorIndexedPath(size_t dep,
                                     const std::vector<CG<Base> >& depVals,
                                     size_t color,
                                     std::set<const OperationNode<Base>*>& indexedOperations) {
            colorIndexedPath(depRef, depVals[dep], color, indexedOperations);
        }

        void uncolor(const std::vector<CG<Base> >& depVals) {
            std::set<size_t>::const_iterator it;
            for (it = dependents.begin(); it != dependents.end(); ++it) {
                uncolor(depVals[*it].getOperationNode());
            }
        }

        static std::set<const OperationNode<Base>*> findOperationsUsingIndependents(OperationNode<Base>& node) {
            std::set<const OperationNode<Base>*> ops;

            clearUsageCount(&node); // reset clear usage count

            findOperationsWithIndeps(node, ops);

            return ops;
        }

        static inline void clearUsageCount(OperationNode<Base>* node) {
            if (node == NULL || node->getUsageCount() == 0)
                return;

            node->resetUsageCount();

            const std::vector<Argument<Base> >& args = node->getArguments();
            size_t size = args.size();
            for (size_t a = 0; a < size; a++) {
                clearUsageCount(args[a].getOperation());
            }
        }

        static inline void uncolor(OperationNode<Base>* node) {
            if (node == NULL || node->getColor() == 0)
                return;

            node->setColor(0);

            const std::vector<Argument<Base> >& args = node->getArguments();
            size_t size = args.size();
            for (size_t a = 0; a < size; a++) {
                uncolor(args[a].getOperation());
            }
        }

        static inline void clearEvaluationOrder(OperationNode<Base>* node) {
            if (node == NULL || node->getEvaluationOrder() == 0)
                return;

            node->setEvaluationOrder(0);

            const std::vector<Argument<Base> >& args = node->getArguments();
            size_t size = args.size();
            for (size_t a = 0; a < size; a++) {
                clearEvaluationOrder(args[a].getOperation());
            }
        }

        inline bool containsConstantIndependent(const OperationNode<Base>* operation, size_t argumentIndex) const {
            typename std::map<const OperationNode<Base>*, std::set<size_t> >::const_iterator it;
            it = constOperationIndependents.find(operation);
            if (it != constOperationIndependents.end()) {
                if (it->second.find(argumentIndex) != it->second.end()) {
                    return true;
                }
            }
            return false;
        }

        /**
         * Determine which independent variables in the loop model do not
         * require an index (always the same for all iterations)
         */
        inline void detectNonIndexedIndependents() {
            typedef typename OperationIndexedIndependents<Base>::MapDep2Indep_type MapIndep2Dep_type;

            // loop operations using independents
            typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::iterator itop2a = indexedOpIndep.op2Arguments.begin();
            while (itop2a != indexedOpIndep.op2Arguments.end()) {
                const OperationNode<Base>* parentOp = itop2a->first;
                OperationIndexedIndependents<Base>& arg2It = itop2a->second;

                bool emptyOp = true;

                // loop the arguments 
                size_t aSize = arg2It.arg2Independents.size();
                for (size_t argIndex = 0; argIndex < aSize; argIndex++) {
                    MapIndep2Dep_type& dep2Ind = arg2It.arg2Independents[argIndex];

                    if (dep2Ind.empty())
                        continue; // argument does not use independents

                    // loop dependents (iterations)
                    bool isIndexed = false;
                    typename MapIndep2Dep_type::const_iterator itDep2Ind = dep2Ind.begin();
                    const OperationNode<Base>* indep = itDep2Ind->second;

                    for (++itDep2Ind; itDep2Ind != dep2Ind.end(); ++itDep2Ind) {
                        if (indep != itDep2Ind->second) {
                            isIndexed = true;
                            break;
                        }
                    }

                    if (!isIndexed) {
                        // make it a non indexed independent
                        constOperationIndependents[parentOp].insert(argIndex);

                        // remove it from the indexed independents
                        dep2Ind.clear();
                    } else {
                        emptyOp = false;
                    }

                }

                if (emptyOp) {
                    indexedOpIndep.op2Arguments.erase(itop2a++);
                } else {
                    ++itop2a;
                }

            }

        }

        virtual ~EquationPattern() {
        }

    private:

        bool comparePath(const CG<Base>& dep1,
                         const CG<Base>& dep2,
                         size_t dep2Index) throw (CGException) {
            if (dep1.getCodeHandler() != dep2.getCodeHandler()) {
                if (dep1.getCodeHandler() != NULL && dep2.getCodeHandler() != NULL)
                    throw CGException("Only one code handler allowed");
                return false;
            }

            if (dep1.isParameter() && dep2.isParameter()) {
                return dep1.getValue() == dep2.getValue();

            } else if (dep1.isVariable() && dep2.isVariable()) {
                OperationNode<Base>* depRefOp = dep1.getOperationNode();
                OperationNode<Base>* dep2Op = dep2.getOperationNode();
                assert(depRefOp->getOperationType() != CGInvOp);

                return comparePath(depRefOp, dep2Op, dep2Index);
            }

            return false;
        }

        bool comparePath(OperationNode<Base>* scRef,
                         OperationNode<Base>* sc2,
                         size_t dep2) {
            saveOperationReference(dep2, sc2, scRef);
            if (dependents.size() == 1) {
                saveOperationReference(depRefIndex, scRef, scRef);
            }

            while (scRef->getOperationType() == CGAliasOp) {
                CPPADCG_ASSERT_KNOWN(scRef->getArguments().size() == 1, "Invalid number of arguments for alias");
                OperationNode<Base>* sc = scRef->getArguments()[0].getOperation();
                if (sc != NULL || sc->getOperationType() == CGInvOp) break;
                scRef = sc;
            }
            while (sc2->getOperationType() == CGAliasOp) {
                CPPADCG_ASSERT_KNOWN(sc2->getArguments().size() == 1, "Invalid number of arguments for alias");
                OperationNode<Base>* sc = sc2->getArguments()[0].getOperation();
                if (sc != NULL || sc->getOperationType() == CGInvOp) break;
                sc2 = sc;
            }

            // check if these nodes where visited before
            if (sc2->getColor() >= minColor_ && scRef->getColor() >= minColor_) {
                /**
                 * been here before for both nodes
                 *  warning: if one would return sc2->getColor() == scRef->getColor()
                 *  it could fail to detect some patterns! e.g.:
                 *    it ref ->  v1 + v1 + v2
                 *    it 2   ->  v3 + v1 + v1
                 *   where v1, v2, v3 have the same expression pattern but 
                 *   correspond to different nodes
                 */
                if (sc2->getColor() == scRef->getColor())
                    return true;
            }
            scRef->setColor(cmpColor_);
            sc2->setColor(cmpColor_);
            cmpColor_++;


            if (scRef->getOperationType() != sc2->getOperationType()) {
                return false;
            }

            assert(scRef->getOperationType() != CGInvOp);

            const std::vector<size_t>& info1 = scRef->getInfo();
            const std::vector<size_t>& info2 = sc2->getInfo();
            if (info1.size() != info2.size()) {
                return false;
            }

            for (size_t e = 0; e < info1.size(); e++) {
                if (info1[e] != info2[e]) {
                    return false;
                }
            }

            const std::vector<Argument<Base> >& args1 = scRef->getArguments();
            const std::vector<Argument<Base> >& args2 = sc2->getArguments();
            size_t size = args1.size();
            if (size != args2.size()) {
                return false;
            }
            for (size_t a = 0; a < size; a++) {
                const Argument<Base>& a1 = args1[a];
                const Argument<Base>& a2 = args2[a];

                if (a1.getParameter() != NULL) {
                    if (a2.getParameter() == NULL || *a1.getParameter() != *a2.getParameter())
                        return false;
                } else {
                    if (a2.getOperation() == NULL) {
                        return false;
                    }
                    OperationNode<Base>* argRefOp = a1.getOperation();
                    OperationNode<Base>* arg2Op = a2.getOperation();
                    bool related;
                    if (argRefOp->getOperationType() == CGInvOp) {
                        related = saveIndependent(scRef, a, argRefOp, arg2Op);
                    } else {
                        related = comparePath(argRefOp, arg2Op, dep2);
                    }

                    if (!related)
                        return false;
                }
            }

            return true; // same pattern
        }

        inline void saveOperationReference(size_t dep2,
                                           const OperationNode<Base>* sc2,
                                           OperationNode<Base>* scRef) {
            operationEO2Reference[dep2][sc2] = scRef;
        }

        bool saveIndependent(const OperationNode<Base>* parentOp,
                             size_t argIndex,
                             const OperationNode<Base>* argRefOp,
                             const OperationNode<Base>* arg2Op) {
            if (argRefOp->getOperationType() != CGInvOp || arg2Op->getOperationType() != CGInvOp) {
                return false;
            }

            /**
             * Must consider that the independent might change from iteration to
             * iteration (even if now it won't)
             */
            typename std::map<const OperationNode<Base>*, std::set<size_t> >::const_iterator it;
            it = constOperationIndependents.find(parentOp);
            if (it != constOperationIndependents.end()) {
                if (it->second.find(argIndex) != it->second.end()) {
                    return false;
                }
            }

            OperationIndexedIndependents<Base>& opIndexedIndep = indexedOpIndep.op2Arguments[parentOp];
            opIndexedIndep.arg2Independents.resize(parentOp != NULL ? parentOp->getArguments().size() : 1);

            std::map<size_t, const OperationNode<Base>*>& dep2Indeps = opIndexedIndep.arg2Independents[argIndex];
            if (dep2Indeps.empty())
                dep2Indeps[depRefIndex] = argRefOp;
            dep2Indeps[currDep_] = arg2Op;

            return true; // same pattern
        }

        inline void colorIndexedPath(const CG<Base>& depRef,
                                     const CG<Base>& dep2,
                                     size_t color,
                                     std::set<const OperationNode<Base>*>& indexedOperations) {
            if (depRef.isVariable() && dep2.isVariable()) {
                OperationNode<Base>* depRefOp = depRef.getOperationNode();
                OperationNode<Base>* dep2Op = dep2.getOperationNode();
                if (depRefOp->getOperationType() != CGInvOp) {
                    colorIndexedPath(depRefOp, dep2Op, color, indexedOperations);
                } else {

                    typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::iterator itop2a;
                    itop2a = indexedOpIndep.op2Arguments.find(NULL);
                    if (itop2a != indexedOpIndep.op2Arguments.end() && !itop2a->second.arg2Independents[0].empty()) {
                        // depends on an index
                        indexedOperations.insert(NULL);
                    }
                }
            }
        }

        inline bool colorIndexedPath(const OperationNode<Base>* scRef,
                                     OperationNode<Base>* sc2,
                                     size_t color,
                                     std::set<const OperationNode<Base>*>& indexedOperations) {

            while (scRef->getOperationType() == CGAliasOp) {
                CPPADCG_ASSERT_KNOWN(scRef->getArguments().size() == 1, "Invalid number of arguments for alias");
                OperationNode<Base>* sc = scRef->getArguments()[0].getOperation();
                if (sc != NULL || sc->getOperationType() == CGInvOp) break;
                scRef = sc;
            }
            while (sc2->getOperationType() == CGAliasOp) {
                CPPADCG_ASSERT_KNOWN(sc2->getArguments().size() == 1, "Invalid number of arguments for alias");
                OperationNode<Base>* sc = sc2->getArguments()[0].getOperation();
                if (sc != NULL || sc->getOperationType() == CGInvOp) break;
                sc2 = sc;
            }

            assert(scRef->getOperationType() == sc2->getOperationType());

            const std::vector<Argument<Base> >& argsRef = scRef->getArguments();

            typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base> >::iterator itop2a;
            bool searched = false;
            bool indexedDependentPath = false;
            bool usesIndexedIndependent = false; // directly uses an indexed independent

            size_t size = argsRef.size();
            for (size_t a = 0; a < size; a++) {
                OperationNode<Base>* argRefOp = argsRef[a].getOperation();
                if (argRefOp != NULL) {
                    bool indexedArg = false;
                    if (argRefOp->getOperationType() == CGInvOp) {
                        // same independent variable can be used in multiple iterations
                        if (!searched) {
                            itop2a = indexedOpIndep.op2Arguments.find(scRef);
                            searched = true;
                        }
                        if (itop2a != indexedOpIndep.op2Arguments.end() && !itop2a->second.arg2Independents[a].empty()) {
                            // depends on an index
                            indexedArg = true;
                            indexedDependentPath = true;
                            usesIndexedIndependent = true;
                        }
                    }

                    if (!indexedArg) {
                        const std::vector<Argument<Base> >& args2 = sc2->getArguments();
                        assert(size == args2.size());
                        indexedDependentPath |= colorIndexedPath(argsRef[a].getOperation(), args2[a].getOperation(), color, indexedOperations);
                    }
                }
            }

            if (indexedDependentPath)
                sc2->setColor(color);
            else
                sc2->setColor(0);

            if (usesIndexedIndependent)
                indexedOperations.insert(sc2);

            return indexedDependentPath;
        }

        static void findOperationsWithIndeps(OperationNode<Base>& node, std::set<const OperationNode<Base>*>& ops) {
            if (node.getUsageCount() > 0)
                return; // been here before

            node.increaseUsageCount();

            const std::vector<Argument<Base> >& args = node.getArguments();
            size_t size = args.size();
            for (size_t a = 0; a < size; a++) {
                OperationNode<Base>* argOp = args[a].getOperation();
                if (argOp != NULL) {
                    if (argOp->getOperationType() == CGInvOp) {
                        ops.insert(&node);
                    } else {
                        findOperationsWithIndeps(*argOp, ops);
                    }
                }
            }
        }

        EquationPattern(const EquationPattern<Base>& other); // not implemented
        EquationPattern& operator=(const EquationPattern<Base>& rhs); // not implemented

    };

}

#endif
