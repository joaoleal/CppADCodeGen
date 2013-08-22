#ifndef CPPAD_CG_EQUATION_PATTERN_INCLUDED
#define CPPAD_CG_EQUATION_PATTERN_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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

namespace CppAD {

    template<class Base>
    class OperationIndexedIndependents {
    public:
        typedef std::map<size_t, const OperationNode<Base>*> MapIndep2Dep_type;
        /**
         * maps the argument index to the several independents used in different
         * equations with the same pattern
         * argument index -> dependent using it -> independent
         */
        std::map<size_t, MapIndep2Dep_type> arg2Independents;
    public:

        inline MapIndep2Dep_type& getDep2Indep(size_t argument) {
            return arg2Independents[argument];
        }
    };

    template<class Base>
    class IndexedIndependent {
    public:
        std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base>* > op2Arguments;
    public:

        OperationIndexedIndependents<Base>* arguments(const OperationNode<Base>* operation) const {
            return op2Arguments.at(operation);
        }

        bool isIndexedOperationArgument(const OperationNode<Base>* node, size_t argIndex) const {
            typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base>*>::const_iterator itIndexes;
            itIndexes = op2Arguments.find(node);
            if (itIndexes == op2Arguments.end()) {
                return false;
            }
            const OperationIndexedIndependents<Base>* indexedArgs = itIndexes->second;
            return indexedArgs->arg2Independents.find(argIndex) != indexedArgs->arg2Independents.end();
        }

        inline OperationIndexedIndependents<Base>* find(const OperationNode<Base>* operation) const {
            typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base>* >::const_iterator it = op2Arguments.find(operation);
            if (it == op2Arguments.end())
                return NULL;
            else
                return it->second;
        }

        //        inline void combine(const IndexedIndependent<Base>& other) {
        //            op2Arguments.insert(other.op2Arguments.begin(), other.op2Arguments.end());
        //        }

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
         * maps operations used by all dependents to the operations of the 
         * reference dependent
         */
        std::map<const OperationNode<Base>*, OperationNode<Base>*> operation2Reference;
        /**
         * Maps the operations that used an indexed independents as direct
         * arguments (reference operation -> argument indexes -> independents)
         */
        IndexedIndependent<Base> indexedOpIndep;
        /**
         * reference operation -> non indexed argument positions
         */
        std::map<const OperationNode<Base>*, std::set<size_t> > constOperationIndependents;

        /**
         * reference independent -> independents
         */
        //std::map<const OperationNode<Base>*, std::vector<const OperationNode<Base>*> > independents;
        //std::set<const OperationNode<Base>*> commonIndependents;
        //std::map<const OperationNode<Base>*, std::set<size_t> > indexedOperationArguments;
    private:
        size_t currDep_;
    public:

        explicit EquationPattern(const CG<Base>& ref,
                                 size_t iDepRef) :
            depRef(ref),
            depRefIndex(iDepRef) {
            dependents.insert(iDepRef);
        }

        bool testAdd(size_t iDep2, const CG<Base>& dep2) {
            IndexedIndependent<Base> independentsBackup = indexedOpIndep;
            std::map<const OperationNode<Base>*, std::set<size_t> > constOperationIndependentsBackup = constOperationIndependents;
            std::map<const OperationNode<Base>*, OperationNode<Base>*> operation2ReferenceBackup = operation2Reference;

            currDep_ = iDep2;

            if (comparePath(depRef, dep2)) {
                dependents.insert(iDep2);
                return true; // matches the reference pattern
            } else {
                // restore
                indexedOpIndep.op2Arguments.swap(independentsBackup.op2Arguments);
                constOperationIndependents.swap(constOperationIndependentsBackup);
                operation2Reference.swap(operation2ReferenceBackup);
                return false; // cannot be added
            }
        }

        inline void colorIndexedPath(size_t dep,
                                     const std::vector<CG<Base> >& depVals,
                                     size_t color,
                                     std::set<const OperationNode<Base>*>& indexedOperations) {
            assert(dependents.size() > 1);
            size_t otherDep;
            if (*dependents.begin() != dep) {
                otherDep = *dependents.begin();
            } else {
                otherDep = *(++dependents.begin());
            }

            colorIndexedPath(depVals[dep], depVals[otherDep], color, indexedOperations);
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

        static inline void uncolor(OperationNode<Base>* sc) {
            if (sc == NULL || sc->getColor() == 0)
                return;

            sc->setColor(0);

            const std::vector<Argument<Base> >& args = sc->getArguments();
            size_t size = args.size();
            for (size_t a = 0; a < size; a++) {
                uncolor(args[a].getOperation());
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

        virtual ~EquationPattern() {
            typename std::map<const OperationNode<Base>*, OperationIndexedIndependents<Base>* >::const_iterator it;

            for (it = indexedOpIndep.op2Arguments.begin(); it != indexedOpIndep.op2Arguments.end(); ++it) {
                delete it->second;
            }

        }

    private:

        bool comparePath(const CG<Base>& dep1,
                         const CG<Base>& dep2) {
            if (dep1.getCodeHandler() != dep2.getCodeHandler()) {
                if (dep1.getCodeHandler() != NULL && dep2.getCodeHandler() != NULL)
                    throw CGException("Only one code handler allowed");
                return false;
            }

            if (dep1.isParameter() && dep2.isParameter()) {
                return dep1.getValue() == dep2.getValue();

            } else if (dep1.isVariable() && dep2.isVariable()) {
                const OperationNode<Base>* depRefOp = dep1.getOperationNode();
                const OperationNode<Base>* dep2Op = dep2.getOperationNode();
                if (depRefOp->getOperationType() == CGInvOp) {
                    return saveIndependent(NULL, 0, depRefOp, dep2Op);
                } else {
                    return comparePath(dep1.getOperationNode(),
                                       dep2.getOperationNode());
                }
            }

            return false;
        }

        bool comparePath(OperationNode<Base>* scRef,
                         OperationNode<Base>* sc2) {
            operation2Reference[sc2] = scRef;
            if (dependents.size() == 1) {
                operation2Reference[scRef] = scRef;
            }

            while (scRef->getOperationType() == CGAliasOp) {
                CPPADCG_ASSERT_KNOWN(scRef->getArguments().size() == 1, "Invalid number of arguments for alias");
                scRef = scRef->getArguments()[0].getOperation();
            }
            while (sc2->getOperationType() == CGAliasOp) {
                CPPADCG_ASSERT_KNOWN(sc2->getArguments().size() == 1, "Invalid number of arguments for alias");
                sc2 = sc2->getArguments()[0].getOperation();
            }

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
                        related = comparePath(argRefOp, arg2Op);
                    }

                    if (!related)
                        return false;
                }
            }

            return true; // same pattern
        }

        bool saveIndependent(const OperationNode<Base>* parentOp,
                             size_t argIndex,
                             const OperationNode<Base>* argRefOp,
                             const OperationNode<Base>* arg2Op) {
            if (argRefOp->getOperationType() != CGInvOp || arg2Op->getOperationType() != CGInvOp) {
                return false;
            }

            if (argRefOp != arg2Op) {
                typename std::map<const OperationNode<Base>*, std::set<size_t> >::const_iterator it;
                it = constOperationIndependents.find(parentOp);
                if (it != constOperationIndependents.end()) {
                    if (it->second.find(argIndex) != it->second.end()) {
                        return false;
                    }
                }

                OperationIndexedIndependents<Base>* opIndexedIndep = indexedOpIndep.op2Arguments[parentOp];
                if (opIndexedIndep == NULL) {
                    opIndexedIndep = new OperationIndexedIndependents<Base>();
                    indexedOpIndep.op2Arguments[parentOp] = opIndexedIndep;
                }

                std::map<size_t, const OperationNode<Base>*>& dep2Indeps = opIndexedIndep->arg2Independents[argIndex];
                if (dep2Indeps.empty())
                    dep2Indeps[depRefIndex] = argRefOp;
                dep2Indeps[currDep_] = arg2Op;
            } else {
                if (indexedOpIndep.isIndexedOperationArgument(parentOp, argIndex)) {
                    return false;
                }
                constOperationIndependents[parentOp].insert(argIndex);
            }

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
                }
            }
        }

        inline bool colorIndexedPath(OperationNode<Base>* scRef,
                                     const OperationNode<Base>* sc2,
                                     size_t color,
                                     std::set<const OperationNode<Base>*>& indexedOperations) {

            while (scRef->getOperationType() == CGAliasOp) {
                CPPADCG_ASSERT_KNOWN(scRef->getArguments().size() == 1, "Invalid number of arguments for alias");
                scRef = scRef->getArguments()[0].getOperation();
            }
            while (sc2->getOperationType() == CGAliasOp) {
                CPPADCG_ASSERT_KNOWN(sc2->getArguments().size() == 1, "Invalid number of arguments for alias");
                sc2 = sc2->getArguments()[0].getOperation();
            }

            if (scRef->getOperationType() == CGInvOp && scRef == sc2) {
                return false; // does NOT depend on an index
            }

            const std::vector<Argument<Base> >& argsRef = scRef->getArguments();
            const std::vector<Argument<Base> >& args2 = sc2->getArguments();
            size_t size = argsRef.size();
            assert(size == args2.size());

            bool indexedDependentPath = false;
            bool usesIndexedIndependent = false; // directly uses an indexed independent
            for (size_t a = 0; a < size; a++) {
                OperationNode<Base>* arg1Op = argsRef[a].getOperation();
                OperationNode<Base>* arg2Op = args2[a].getOperation();
                if (arg1Op != NULL) {
                    if (arg1Op->getOperationType() == CGInvOp && arg1Op != arg2Op) {
                        // depends on an index
                        indexedDependentPath = true;
                        usesIndexedIndependent = true;
                    } else {
                        indexedDependentPath |= colorIndexedPath(argsRef[a].getOperation(), args2[a].getOperation(), color, indexedOperations);
                    }
                }
            }

            if (indexedDependentPath)
                scRef->setColor(color);
            else
                scRef->setColor(0);

            if (usesIndexedIndependent)
                indexedOperations.insert(scRef);

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
