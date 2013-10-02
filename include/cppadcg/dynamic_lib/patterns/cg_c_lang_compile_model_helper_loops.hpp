#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_INCLUDED
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

    namespace loops {

        /***********************************************************************
         *  Utility classes
         **********************************************************************/

        template<class Base>
        class IfBranchInfo {
        public:
            std::set<size_t> iterations;
            OperationNode<Base>* node;
        };

        template <class Base>
        class IfElseInfo {
        public:
            std::map<SizeN1stIt, IfBranchInfo<Base> > firstIt2Branch;
            OperationNode<Base>* endIf;

            inline IfElseInfo() :
                endIf(NULL) {
            }
        };

        template<class Base>
        class IndexedDependentLoopInfo {
        public:
            std::vector<size_t> indexes;
            std::vector<CG<Base> > origVals;
            IndexPattern* pattern;

            inline IndexedDependentLoopInfo() :
                pattern(NULL) {
            }
        };

        template<class Base>
        inline vector<CG<Base> > createIndexedIndependents(CodeHandler<Base>& handler,
                                                           LoopModel<Base>& loop,
                                                           IndexOperationNode<Base>& iterationIndexOp) {

            const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = loop.getIndexedIndepIndexes();
            size_t nIndexed = indexedIndepIndexes.size();

            vector<CG<Base> > x(nIndexed); // zero order

            std::vector<Argument<Base> > xIndexedArgs(1);
            xIndexedArgs[0] = Argument<Base>(iterationIndexOp);
            std::vector<size_t> info(2);
            info[0] = 0; // tx

            for (size_t j = 0; j < nIndexed; j++) {
                info[1] = handler.addLoopIndependentIndexPattern(*loop.getIndependentIndexPatterns()[j], j);
                x[j] = handler.createCG(new OperationNode<Base>(CGLoopIndexedIndepOp, info, xIndexedArgs));
            }

            return x;
        }

        template<class Base>
        inline vector<CG<Base> > createLoopIndependentVector(CodeHandler<Base>& handler,
                                                             LoopModel<Base>& loop,
                                                             const vector<CG<Base> >& indexedIndeps,
                                                             const vector<CG<Base> >& nonIndexed,
                                                             const vector<CG<Base> >& nonIndexedTmps) {

            const std::vector<std::vector<LoopPosition> >& indexedIndepIndexes = loop.getIndexedIndepIndexes();
            const std::vector<LoopPosition>& nonIndexedIndepIndexes = loop.getNonIndexedIndepIndexes();
            const std::vector<LoopPosition>& temporaryIndependents = loop.getTemporaryIndependents();

            size_t nIndexed = indexedIndepIndexes.size();
            size_t nNonIndexed = nonIndexedIndepIndexes.size();
            size_t nTape = nIndexed + nNonIndexed + temporaryIndependents.size();

            // indexed independents
            vector<CG<Base> > x(nTape);
            for (size_t j = 0; j < nIndexed; j++) {
                x[j] = indexedIndeps[j];
            }

            // non indexed
            for (size_t j = 0; j < nNonIndexed; j++) {
                x[nIndexed + j] = nonIndexed[nonIndexedIndepIndexes[j].original];
            }

            // temporaries
            for (size_t j = 0; j < temporaryIndependents.size(); j++) {
                x[nIndexed + nNonIndexed + j] = nonIndexedTmps[temporaryIndependents[j].original];
            }

            return x;
        }

        template<class Base>
        inline vector<CG<Base> > createLoopDependentVector(CodeHandler<Base>& handler,
                                                           LoopModel<Base>& loop,
                                                           IndexOperationNode<Base>& iterationIndexOp) {

            const vector<IndexPattern*>& depIndexes = loop.getDependentIndexPatterns();
            vector<CG<Base> > deps(depIndexes.size());

            size_t dep_size = depIndexes.size();
            size_t x_size = loop.getTapeIndependentCount();

            std::vector<Argument<Base> > xIndexedArgs(1);
            xIndexedArgs[0] = Argument<Base>(iterationIndexOp);
            std::vector<size_t> info(2);
            info[0] = 2; // py2

            for (size_t i = 0; i < dep_size; i++) {
                IndexPattern* ip = depIndexes[i];
                info[1] = handler.addLoopIndependentIndexPattern(*ip, x_size + i); // dependent index pattern location
                deps[i] = handler.createCG(new OperationNode<Base>(CGLoopIndexedIndepOp, info, xIndexedArgs));
            }

            return deps;
        }

        /***************************************************************************
         *  Methods related with loop insertion into the operation graph
         **************************************************************************/

        template<class Base>
        LoopEndOperationNode<Base>* createLoopEnd(CodeHandler<Base>& handler,
                                                  LoopStartOperationNode<Base>& loopStart,
                                                  const vector<std::pair<CG<Base>, IndexPattern*> >& indexedLoopResults,
                                                  const std::set<IndexOperationNode<Base>*>& indexesOps,
                                                  size_t assignOrAdd) {
            std::vector<Argument<Base> > endArgs;
            std::vector<Argument<Base> > indexedArgs(1 + indexesOps.size());
            std::vector<size_t> info(2);

            size_t dep_size = indexedLoopResults.size();
            endArgs.reserve(dep_size);

            for (size_t i = 0; i < dep_size; i++) {
                const std::pair<CG<Base>, IndexPattern*>& depInfo = indexedLoopResults[i];
                if (depInfo.second != NULL) {
                    indexedArgs.resize(1);

                    indexedArgs[0] = asArgument(depInfo.first); // indexed expression
                    typename std::set<IndexOperationNode<Base>*>::const_iterator itIndexOp;
                    for (itIndexOp = indexesOps.begin(); itIndexOp != indexesOps.end(); ++itIndexOp) {
                        indexedArgs.push_back(Argument<Base>(**itIndexOp)); // dependency on the index
                    }

                    info[0] = handler.addLoopDependentIndexPattern(*depInfo.second); // dependent index pattern location
                    info[1] = assignOrAdd;

                    OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, info, indexedArgs);
                    handler.manageOperationNodeMemory(yIndexed);
                    endArgs.push_back(Argument<Base>(*yIndexed));
                } else {
                    OperationNode<Base>* n = depInfo.first.getOperationNode();
                    assert(n != NULL);
                    endArgs.push_back(Argument<Base>(*n));
                }
            }

            LoopEndOperationNode<Base>* loopEnd = new LoopEndOperationNode<Base>(loopStart, endArgs);
            handler.manageOperationNodeMemory(loopEnd);

            return loopEnd;
        }

        template<class Base>
        inline void moveNonIndexedOutsideLoop(LoopStartOperationNode<Base>& loopStart,
                                              LoopEndOperationNode<Base>& loopEnd) {
            //EquationPattern<Base>::uncolor(dependents[dep].getOperationNode());
            const IndexDclrOperationNode<Base>& loopIndex = loopStart.getIndex();
            std::set<OperationNode<Base>*> nonIndexed;

            const std::vector<Argument<Base> >& endArgs = loopEnd.getArguments();
            for (size_t i = 0; i < endArgs.size(); i++) {
                assert(endArgs[i].getOperation() != NULL);
                findNonIndexedNodes(*endArgs[i].getOperation(), nonIndexed, loopIndex);
            }

            std::vector<Argument<Base> >& startArgs = loopStart.getArguments();

            size_t sas = startArgs.size();
            startArgs.resize(sas + nonIndexed.size());
            size_t i = 0;
            typename std::set<OperationNode<Base>*>::const_iterator it;
            for (it = nonIndexed.begin(); it != nonIndexed.end(); ++it, i++) {
                startArgs[sas + i] = Argument<Base>(**it);
            }
        }

        template<class Base>
        inline bool findNonIndexedNodes(OperationNode<Base>& node,
                                        std::set<OperationNode<Base>*>& nonIndexed,
                                        const IndexDclrOperationNode<Base>& loopIndex) {
            if (node.getColor() > 0)
                return node.getColor() == 1;

            if (node.getOperationType() == CGIndexDeclarationOp) {
                if (&node == &loopIndex) {
                    node.setColor(2);
                    return false; // depends on the loop index
                }
            }

            const std::vector<Argument<Base> >& args = node.getArguments();
            size_t size = args.size();

            bool indexedPath = false; // whether or not this node depends on indexed independents
            bool nonIndexedArgs = false; // whether or not there are non indexed arguments
            for (size_t a = 0; a < size; a++) {
                OperationNode<Base>* arg = args[a].getOperation();
                if (arg != NULL) {
                    bool nonIndexedArg = findNonIndexedNodes(*arg, nonIndexed, loopIndex);
                    nonIndexedArgs |= nonIndexedArg;
                    indexedPath |= !nonIndexedArg;
                }
            }

            node.setColor(indexedPath ? 2 : 1);

            if (node.getOperationType() == CGArrayElementOp ||
                    node.getOperationType() == CGAtomicForwardOp ||
                    node.getOperationType() == CGAtomicReverseOp) {
                return !indexedPath; // should not move array creation elements outside the loop
            }

            if (indexedPath && nonIndexedArgs) {
                for (size_t a = 0; a < size; a++) {
                    OperationNode<Base>* arg = args[a].getOperation();
                    if (arg != NULL && arg->getColor() == 1 && // must be a non indexed expression
                            arg->getOperationType() != CGInvOp) {// no point in moving just one variable outside

                        nonIndexed.insert(arg);
                    }
                }
            }

            return !indexedPath;
        }

        template<class Base>
        inline IfElseInfo<Base>* findExistingIfElse(vector<IfElseInfo<Base> >& ifElses,
                                                    const std::map<SizeN1stIt, std::pair<size_t, std::set<size_t> > >& first2Iterations) {
            using namespace std;

            // try to find an existing if-else where these operations can be added
            for (size_t f = 0; f < ifElses.size(); f++) {
                IfElseInfo<Base>& ifElse = ifElses[f];

                if (first2Iterations.size() != ifElse.firstIt2Branch.size())
                    continue;

                bool matches = true;
                map<SizeN1stIt, pair<size_t, set<size_t> > >::const_iterator itLoc = first2Iterations.begin();
                typename map<SizeN1stIt, IfBranchInfo<Base> >::const_iterator itBranches = ifElse.firstIt2Branch.begin();
                for (; itLoc != first2Iterations.end(); ++itLoc, ++itBranches) {
                    if (itLoc->second.second != itBranches->second.iterations) {
                        matches = false;
                        break;
                    }
                }

                if (matches) {
                    return &ifElse;
                }
            }

            return NULL;
        }
    }
}

#endif
