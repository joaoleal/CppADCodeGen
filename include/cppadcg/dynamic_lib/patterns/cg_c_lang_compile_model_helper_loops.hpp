#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_INCLUDED
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

    namespace loops {

        /***********************************************************************
         *  Utility classes
         **********************************************************************/
        template<class Base>
        class IfBranchData {
        public:
            CG<Base> value;
            std::map<size_t, size_t> locations;
        public:

            inline IfBranchData() {
            }

            inline IfBranchData(const CG<Base>& v,
                                const std::map<size_t, size_t>& loc) :
                value(v),
                locations(loc) {
            }
        };

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

        class JacobianWithLoopsRowInfo {
        public:
            // tape J index -> {locationIt0, locationIt1, ...}
            std::map<size_t, std::vector<size_t> > indexedPositions;

            // original J index -> {locationIt0, locationIt1, ...}
            std::map<size_t, std::vector<size_t> > nonIndexedPositions;

            // original J index 
            std::set<size_t> nonIndexedEvals;

            // original J index -> k index 
            std::map<size_t, std::set<size_t> > tmpEvals;
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

        template<class Base>
        inline CG<Base> createLoopDependentFunctionResult(CodeHandler<Base>& handler,
                                                          size_t i, const CG<Base>& val, IndexPattern* ip,
                                                          IndexOperationNode<Base>& iterationIndexOp) {

            size_t assignOrAdd = 1; // add
            std::vector<Argument<Base> > indexedArgs(2);
            std::vector<size_t> aInfo(2);

            if (ip != NULL) {
                aInfo[0] = handler.addLoopDependentIndexPattern(*ip); // dependent index pattern location
                aInfo[1] = assignOrAdd;
                indexedArgs[0] = asArgument(val); // indexed expression
                indexedArgs[1] = Argument<Base>(iterationIndexOp); // index  ///jrowIndexOp

                OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, aInfo, indexedArgs);
                handler.manageOperationNodeMemory(yIndexed);

                return handler.createCG(Argument<Base>(*yIndexed));

            } else if (val.getOperationNode() != NULL &&
                    val.getOperationNode()->getOperationType() == CGEndIfOp) {

                std::vector<size_t> info(1);
                info[0] = i; // points to itself
                std::vector<Argument<Base> > args(1);
                args[0] = Argument<Base>(*val.getOperationNode());

                return handler.createCG(new OperationNode<Base> (CGDependentRefRhsOp, info, args));

            } else {
                return val;
            }
        }

        /***************************************************************************
         *  Methods related with loop insertion into the operation graph
         **************************************************************************/

        template<class Base>
        LoopEndOperationNode<Base>* createLoopEnd(CodeHandler<Base>& handler,
                                                  LoopStartOperationNode<Base>& loopStart,
                                                  const std::vector<std::pair<CG<Base>, IndexPattern*> >& indexedLoopResults,
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
                    CPPADCG_ASSERT_UNKNOWN(n != NULL);
                    endArgs.push_back(Argument<Base>(*n));
                }
            }

            LoopEndOperationNode<Base>* loopEnd = new LoopEndOperationNode<Base>(loopStart, endArgs);
            handler.manageOperationNodeMemory(loopEnd);

            return loopEnd;
        }

        template<class Base>
        inline void moveNonIndexedOutsideLoop(CodeHandler<Base>& handler,
                                              LoopStartOperationNode<Base>& loopStart,
                                              LoopEndOperationNode<Base>& loopEnd) {
            //EquationPattern<Base>::uncolor(dependents[dep].getOperationNode());
            const IndexDclrOperationNode<Base>& loopIndex = loopStart.getIndex();
            std::set<OperationNode<Base>*> nonIndexed;

            const std::vector<Argument<Base> >& endArgs = loopEnd.getArguments();
            for (size_t i = 0; i < endArgs.size(); i++) {
                CPPADCG_ASSERT_UNKNOWN(endArgs[i].getOperation() != NULL);
                LoopNonIndexedLocator<Base>(handler, nonIndexed, loopIndex).findNonIndexedNodes(*endArgs[i].getOperation());
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
        class LoopNonIndexedLocator {
        private:
            CodeHandler<Base>& handler_;
            std::set<OperationNode<Base>*>& nonIndexed_;
            const IndexDclrOperationNode<Base>& loopIndex_;
        public:

            inline LoopNonIndexedLocator(CodeHandler<Base>& handler,
                                         std::set<OperationNode<Base>*>& nonIndexed,
                                         const IndexDclrOperationNode<Base>& loopIndex) :
                handler_(handler),
                nonIndexed_(nonIndexed),
                loopIndex_(loopIndex) {
            }

            inline bool findNonIndexedNodes(OperationNode<Base>& node) {
                if (node.getColor() > 0)
                    return node.getColor() == 1;

                if (node.getOperationType() == CGIndexDeclarationOp) {
                    if (&node == &loopIndex_) {
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
                        bool nonIndexedArg = findNonIndexedNodes(*arg);
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
                        if (arg != NULL && arg->getColor() == 1) {// must be a non indexed expression
                            CGOpCode op = arg->getOperationType();
                            if (op != CGInvOp && op != CGTmpDclOp) {// no point in moving just one variable outside

                                if (op == CGLoopIndexedTmpOp) {
                                    // must not place a CGLoopIndexedTmpOp operation outside the loop
                                    Argument<Base> assignArg = arg->getArguments()[1];
                                    if (assignArg.getOperation() != NULL) { // no point in moving a constant value outside
                                        OperationNode<Base>* assignNode = new OperationNode<Base>(CGAssignOp, assignArg);
                                        handler_.manageOperationNodeMemory(assignNode);
                                        arg->getArguments()[1] = Argument<Base>(*assignNode);
                                        nonIndexed_.insert(assignNode);
                                    }
                                } else {
                                    nonIndexed_.insert(arg);
                                }
                            }
                        }
                    }
                }

                return !indexedPath;
            }
        };

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

        template<class Base>
        OperationNode<Base>* createIndexConditionExpressionOp(const std::set<size_t>& iterations,
                                                              const std::set<size_t>& usedIter,
                                                              size_t maxIter,
                                                              IndexOperationNode<Base>& iterationIndexOp) {
            std::vector<size_t> info = createIndexConditionExpression(iterations, usedIter, maxIter);
            std::vector<Argument<Base> > args(1);
            args[0] = Argument<Base>(iterationIndexOp);

            return new OperationNode<Base>(CGIndexCondExprOp, info, args);
        }

        std::vector<size_t> createIndexConditionExpression(const std::set<size_t>& iterations,
                                                           const std::set<size_t>& usedIter,
                                                           size_t maxIter) {
            CPPADCG_ASSERT_UNKNOWN(!iterations.empty());

            std::map<size_t, bool> allIters;
            for (std::set<size_t>::const_iterator it = usedIter.begin(); it != usedIter.end(); ++it) {
                allIters[*it] = false;
            }
            for (std::set<size_t>::const_iterator it = iterations.begin(); it != iterations.end(); ++it) {
                allIters[*it] = true;
            }

            std::vector<size_t> info;
            info.reserve(iterations.size() / 2 + 2);

            std::map<size_t, bool>::const_iterator it = allIters.begin();
            while (it != allIters.end()) {
                std::map<size_t, bool> ::const_iterator min = it;
                std::map<size_t, bool>::const_iterator max = it;
                std::map<size_t, bool>::const_iterator minNew = allIters.end();
                std::map<size_t, bool>::const_iterator maxNew = allIters.end();
                if (it->second) {
                    minNew = it;
                    maxNew = it;
                }

                for (++it; it != allIters.end(); ++it) {
                    if (it->first != max->first + 1) {
                        break;
                    }

                    max = it;
                    if (it->second) {
                        if (minNew == allIters.end())
                            minNew = it;
                        maxNew = it;
                    }
                }

                if (minNew != allIters.end()) {
                    // contains elements from the current iteration set
                    if (maxNew->first == minNew->first) {
                        // only one element
                        info.push_back(minNew->first);
                        info.push_back(maxNew->first);
                    } else {
                        //several elements
                        if (min->first == 0)
                            info.push_back(min->first);
                        else
                            info.push_back(minNew->first);

                        if (max->first == maxIter)
                            info.push_back(std::numeric_limits<size_t>::max());
                        else
                            info.push_back(maxNew->first);
                    }
                }
            }

            return info;
        }

        template<class Base>
        inline CG<Base> createConditionalContribution(CodeHandler<Base>& handler,
                                                      const std::map<size_t, IfBranchData<Base> >& branches,
                                                      size_t maxIter,
                                                      size_t nLocalIter,
                                                      IndexOperationNode<Base>& iterationIndexOp,
                                                      vector<IfElseInfo<Base> >& ifElses,
                                                      bool printResult = false) {
            using namespace std;

            map<SizeN1stIt, pair<size_t, set<size_t> > > firstIt2Count2Iterations;
            typename map<size_t, IfBranchData<Base> >::const_iterator itb;
            for (itb = branches.begin(); itb != branches.end(); ++itb) {
                set<size_t> iterations;
                mapKeys(itb->second.locations, iterations);
                firstIt2Count2Iterations[SizeN1stIt(iterations.size(), *iterations.begin())] = make_pair(itb->first, iterations);
            }

            IfElseInfo<Base>* ifElseBranches = findExistingIfElse(ifElses, firstIt2Count2Iterations);
            bool reusingIfElse = ifElseBranches != NULL;
            if (!reusingIfElse) {
                size_t s = ifElses.size();
                ifElses.resize(s + 1);
                ifElseBranches = &ifElses[s];
            }

            /**
             * create/change each if/else branch
             */
            OperationNode<Base>* ifStart = NULL;
            OperationNode<Base>* ifBranch = NULL;
            Argument<Base> nextBranchArg;
            set<size_t> usedIter;

            map<SizeN1stIt, pair<size_t, set<size_t> > >::const_iterator it1st2Count2Iters;
            for (it1st2Count2Iters = firstIt2Count2Iterations.begin(); it1st2Count2Iters != firstIt2Count2Iterations.end(); ++it1st2Count2Iters) {
                size_t firstIt = it1st2Count2Iters->first.second;
                size_t count = it1st2Count2Iters->second.first;
                const set<size_t>& iterations = it1st2Count2Iters->second.second;
                const IfBranchData<Base>& branchData = branches.at(count);

                size_t iterCount = iterations.size();

                SizeN1stIt pos(iterCount, firstIt);

                if (reusingIfElse) {
                    //reuse existing node
                    ifBranch = ifElseBranches->firstIt2Branch.at(pos).node;
                    if (nextBranchArg.getOperation() != NULL)
                        ifBranch->getArguments().push_back(nextBranchArg);

                } else if (usedIter.size() + iterCount == nLocalIter) {
                    // all other iterations: ELSE
                    ifBranch = new OperationNode<Base>(CGElseOp, Argument<Base>(*ifBranch), nextBranchArg);
                    handler.manageOperationNodeMemory(ifBranch);
                } else {
                    // depends on the iteration index
                    OperationNode<Base>* cond = createIndexConditionExpressionOp<Base>(iterations, usedIter, maxIter, iterationIndexOp);
                    handler.manageOperationNodeMemory(cond);

                    if (ifStart == NULL) {
                        // IF
                        ifStart = new OperationNode<Base>(CGStartIfOp, Argument<Base>(*cond));
                        ifBranch = ifStart;
                    } else {
                        // ELSE IF
                        ifBranch = new OperationNode<Base>(CGElseIfOp, Argument<Base>(*ifBranch), Argument<Base>(*cond), nextBranchArg);
                    }

                    handler.manageOperationNodeMemory(ifBranch);

                    usedIter.insert(iterations.begin(), iterations.end());
                }

                IndexPattern* pattern = IndexPattern::detect(branchData.locations);
                handler.manageLoopDependentIndexPattern(pattern);

                Argument<Base> value;
                if (printResult) {
                    PrintOperationNode<Base>* printNode = new PrintOperationNode<Base>("__________", asArgument(branchData.value), "\n");
                    handler.manageOperationNodeMemory(printNode);
                    value = Argument<Base>(*printNode);
                } else {
                    value = asArgument(branchData.value);
                }

                std::vector<size_t> ainfo(2);
                ainfo[0] = handler.addLoopDependentIndexPattern(*pattern); // dependent index pattern location
                ainfo[1] = 1; // assignOrAdd
                std::vector<Argument<Base> > indexedArgs(2);
                indexedArgs[0] = value; // indexed expression
                indexedArgs[1] = Argument<Base>(iterationIndexOp); // dependency on the index
                OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, ainfo, indexedArgs);
                handler.manageOperationNodeMemory(yIndexed);

                OperationNode<Base>* ifAssign = new OperationNode<Base>(CGCondResultOp, Argument<Base>(*ifBranch), Argument<Base>(*yIndexed));
                handler.manageOperationNodeMemory(ifAssign);
                nextBranchArg = Argument<Base>(*ifAssign);

                if (!reusingIfElse) {
                    IfBranchInfo<Base>& branch = ifElseBranches->firstIt2Branch[pos]; // creates a new if branch
                    branch.iterations = iterations;
                    branch.node = ifBranch;
                }
            }

            /**
             * end if
             */
            if (reusingIfElse) {
                ifElseBranches->endIf->getArguments().push_back(nextBranchArg);
            } else {
                ifElseBranches->endIf = new OperationNode<Base>(CGEndIfOp, Argument<Base>(*ifBranch), nextBranchArg);
                handler.manageOperationNodeMemory(ifElseBranches->endIf);
            }

            return handler.createCG(Argument<Base>(*ifElseBranches->endIf));
        }

        /**
         * Contribution to a constant location
         */
        template<class Base>
        CG<Base> createConditionalContribution(CodeHandler<Base>& handler,
                                               LinearIndexPattern& pattern,
                                               const std::set<size_t>& iterations,
                                               size_t maxIter,
                                               const CG<Base>& ddfdxdx,
                                               IndexOperationNode<Base>& iterationIndexOp,
                                               vector<IfElseInfo<Base> >& ifElses) {
            using namespace std;

            CPPADCG_ASSERT_UNKNOWN(pattern.getLinearSlopeDy() == 0); // must be a constant index

            // try to find an existing if-else where these operations can be added
            map<SizeN1stIt, pair<size_t, set<size_t> > > firstIt2Count2Iterations;
            SizeN1stIt pos(iterations.size(), *iterations.begin());
            firstIt2Count2Iterations[pos] = make_pair(1, iterations);

            IfElseInfo<Base>* ifElseBranches = findExistingIfElse(ifElses, firstIt2Count2Iterations);
            bool reusingIfElse = ifElseBranches != NULL;
            if (!reusingIfElse) {
                size_t s = ifElses.size();
                ifElses.resize(s + 1);
                ifElseBranches = &ifElses[s];
            }

            /**
             * create/change each if/else branch
             */
            OperationNode<Base>* ifBranch = NULL;

            if (reusingIfElse) {
                //reuse existing node
                ifBranch = ifElseBranches->firstIt2Branch.at(pos).node;

            } else {
                // depends on the iterations indexes
                const set<size_t> usedIter;
                OperationNode<Base>* cond = createIndexConditionExpressionOp<Base>(iterations, usedIter, maxIter, iterationIndexOp);
                handler.manageOperationNodeMemory(cond);

                ifBranch = new OperationNode<Base>(CGStartIfOp, Argument<Base>(*cond));
                handler.manageOperationNodeMemory(ifBranch);
            }

            std::vector<size_t> ainfo(2);
            ainfo[0] = handler.addLoopDependentIndexPattern(pattern); // dependent index pattern location
            ainfo[1] = 1; // assignOrAdd
            std::vector<Argument<Base> > indexedArgs(2);
            indexedArgs[0] = asArgument(ddfdxdx); // indexed expression
            indexedArgs[1] = Argument<Base>(iterationIndexOp); // dependency on the index
            OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, ainfo, indexedArgs);
            handler.manageOperationNodeMemory(yIndexed);

            OperationNode<Base>* ifAssign = new OperationNode<Base>(CGCondResultOp, Argument<Base>(*ifBranch), Argument<Base>(*yIndexed));
            handler.manageOperationNodeMemory(ifAssign);
            Argument<Base> nextBranchArg = Argument<Base>(*ifAssign);

            if (!reusingIfElse) {
                IfBranchInfo<Base>& branch = ifElseBranches->firstIt2Branch[pos]; // creates a new if branch
                branch.iterations = iterations;
                branch.node = ifBranch;
            }

            /**
             * end if
             */
            if (reusingIfElse) {
                ifElseBranches->endIf->getArguments().push_back(nextBranchArg);
            } else {
                ifElseBranches->endIf = new OperationNode<Base>(CGEndIfOp, Argument<Base>(*ifBranch), nextBranchArg);
                handler.manageOperationNodeMemory(ifElseBranches->endIf);
            }

            return handler.createCG(Argument<Base>(*ifElseBranches->endIf));
        }

        /**
         * 
         * @param handler source code handler
         * @param locationsIter2Pos maps each iteration to the location of the result
         * @param iterCount the number of iteration of the loop
         * @param value the value determined inside the loop
         * @param pattern the pattern used to save the value
         * @param assignOrAdd whether the value is assigned or added to the
         *                    dependent array
         * @param iterationIndexOp the iteration index operation for this loop
         * @param ifElses conditions used inside the loop
         * @return 
         */
        template<class Base>
        std::pair<CG<Base>, IndexPattern*> createLoopResult(CodeHandler<Base>& handler,
                                                            const std::map<size_t, size_t>& locationsIter2Pos,
                                                            size_t iterCount,
                                                            const CG<Base>& value,
                                                            IndexPattern* pattern,
                                                            size_t assignOrAdd,
                                                            IndexOperationNode<Base>& iterationIndexOp,
                                                            vector<IfElseInfo<Base> >& ifElses) {
            using namespace std;

            if (locationsIter2Pos.size() == iterCount) {
                // present in all iterations

                return make_pair(value, pattern);

            } else {
                /**
                 * must create a conditional element so that this 
                 * contribution is only evaluated at the relevant iterations
                 */

                // try to find an existing if-else where these operations can be added
                map<SizeN1stIt, pair<size_t, set<size_t> > > firstIt2Count2Iterations;

                set<size_t> iterations;
                mapKeys(locationsIter2Pos, iterations);
                SizeN1stIt pos(iterations.size(), *iterations.begin());
                firstIt2Count2Iterations[pos] = make_pair(*iterations.begin(), iterations);

                IfElseInfo<Base>* ifElseBranches = findExistingIfElse(ifElses, firstIt2Count2Iterations);
                bool reusingIfElse = ifElseBranches != NULL;
                if (!reusingIfElse) {
                    size_t s = ifElses.size();
                    ifElses.resize(s + 1);
                    ifElseBranches = &ifElses[s];
                }

                OperationNode<Base>* ifStart;

                if (reusingIfElse) {
                    //reuse existing node
                    ifStart = ifElseBranches->firstIt2Branch.at(pos).node;
                } else {
                    // depends on the iterations indexes
                    set<size_t> usedIter;
                    OperationNode<Base>* cond = createIndexConditionExpressionOp<Base>(iterations, usedIter, iterCount - 1, iterationIndexOp);
                    handler.manageOperationNodeMemory(cond);

                    ifStart = new OperationNode<Base>(CGStartIfOp, Argument<Base>(*cond));
                    handler.manageOperationNodeMemory(ifStart);
                }

                std::vector<size_t> ainfo(2);
                ainfo[0] = handler.addLoopDependentIndexPattern(*pattern); // dependent index pattern location
                ainfo[1] = assignOrAdd;
                std::vector<Argument<Base> > indexedArgs(2);
                indexedArgs[0] = asArgument(value); // indexed expression
                indexedArgs[1] = Argument<Base>(iterationIndexOp); // dependency on the index
                OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, ainfo, indexedArgs);
                handler.manageOperationNodeMemory(yIndexed);

                OperationNode<Base>* ifAssign = new OperationNode<Base>(CGCondResultOp, Argument<Base>(*ifStart), Argument<Base>(*yIndexed));
                handler.manageOperationNodeMemory(ifAssign);

                if (!reusingIfElse) {
                    // existing 'if' with the same iterations
                    IfBranchInfo<Base>& branch = ifElseBranches->firstIt2Branch[pos]; // creates a new if branch
                    branch.iterations = iterations;
                    branch.node = ifStart;
                }

                if (reusingIfElse) {
                    ifElseBranches->endIf->getArguments().push_back(Argument<Base>(*ifAssign));
                } else {
                    ifElseBranches->endIf = new OperationNode<Base>(CGEndIfOp, Argument<Base>(*ifStart), Argument<Base>(*ifAssign));
                    handler.manageOperationNodeMemory(ifElseBranches->endIf);
                }

                IndexPattern* p = NULL;
                return make_pair(handler.createCG(Argument<Base>(*ifElseBranches->endIf)), p);
            }

        }

        class ArrayElementCopyPattern {
        public:
            IndexPattern* resultPattern;
            IndexPattern* compressedPattern;
        public:

            inline ArrayElementCopyPattern() :
                resultPattern(NULL),
                compressedPattern(NULL) {
            }

            inline ArrayElementCopyPattern(IndexPattern* resultPat,
                                           IndexPattern* compressedPat) :
                resultPattern(resultPat),
                compressedPattern(compressedPat) {
            }

            inline ~ArrayElementCopyPattern() {
                delete resultPattern;
                delete compressedPattern;
            }

        };

        class ArrayElementGroup {
        public:
            std::set<size_t> keys;
            std::vector<ArrayElementCopyPattern> elements;

            ArrayElementGroup(const std::set<size_t>& k, size_t size) :
                keys(k),
                elements(size) {
            }
        };

        class ArrayGroup {
        public:
            std::auto_ptr<IndexPattern> pattern;
            std::auto_ptr<IndexPattern> startLocPattern;
            SmartMapValuePointer<size_t, ArrayElementGroup> elCount2elements;
        };

        /**
         * 
         * @param loopGroups Used elements from the arrays provided by the group
         *                   function calls (loop->group->{array->{compressed position} })
         * @param userElLocation maps each element to its position 
         *                      (array -> [compressed elements { original index }] )
         * @param orderedArray
         * @param loopCalls
         * @param garbage holds created ArrayGroups
         */
        template<class Base>
        inline void determineForRevUsagePatterns(const std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >& loopGroups,
                                                 const std::map<size_t, std::vector<std::set<size_t> > >& userElLocation,
                                                 const std::map<size_t, bool>& orderedArray,
                                                 std::map<size_t, std::map<LoopModel<Base>*, std::map<size_t, ArrayGroup*> > >& loopCalls,
                                                 SmartVectorPointer<ArrayGroup>& garbage) {

            using namespace std;

            std::vector<size_t> arrayStart;
            /**
             * Determine jcol index patterns and
             * array start patterns
             */
            std::vector<size_t> localit2jcols;
            typename map<LoopModel<Base>*, map<size_t, map<size_t, set<size_t> > > >::const_iterator itlge;
            for (itlge = loopGroups.begin(); itlge != loopGroups.end(); ++itlge) {
                LoopModel<Base>* loop = itlge->first;

                garbage.v.reserve(garbage.v.size() + itlge->second.size());

                map<size_t, map<size_t, set<size_t> > >::const_iterator itg;
                for (itg = itlge->second.begin(); itg != itlge->second.end(); ++itg) {
                    size_t group = itg->first;
                    const map<size_t, set<size_t> >& jcols2e = itg->second;

                    // group by number of iterations
                    std::auto_ptr<ArrayGroup> data(new ArrayGroup());

                    /**
                     * jcol pattern
                     */
                    mapKeys(jcols2e, localit2jcols);
                    data->pattern.reset(IndexPattern::detect(localit2jcols));

                    /**
                     * array start pattern
                     */
                    bool ordered = true;
                    for (size_t l = 0; l < localit2jcols.size(); l++) {
                        if (!orderedArray.at(localit2jcols[l])) {
                            ordered = false;
                            break;
                        }
                    }

                    if (ordered) {
                        arrayStart.resize(localit2jcols.size());

                        for (size_t l = 0; l < localit2jcols.size(); l++) {
                            const std::vector<set<size_t> >& location = userElLocation.at(localit2jcols[l]);
                            arrayStart[l] = *location[0].begin();
                        }

                        data->startLocPattern.reset(IndexPattern::detect(arrayStart));
                    } else {
                        /**
                         * combine calls to this group function which provide 
                         * the same number of elements
                         */
                        map<size_t, map<size_t, size_t> > elCount2localIt2jcols;

                        map<size_t, set<size_t> >::const_iterator itJcols2e;
                        size_t localIt = 0;
                        for (itJcols2e = jcols2e.begin(); itJcols2e != jcols2e.end(); ++itJcols2e, localIt++) {
                            size_t elCount = itJcols2e->second.size();
                            elCount2localIt2jcols[elCount][localIt] = itJcols2e->first;
                        }

                        map<size_t, map<size_t, size_t> > ::const_iterator elC2jcolIt;
                        for (elC2jcolIt = elCount2localIt2jcols.begin(); elC2jcolIt != elCount2localIt2jcols.end(); ++elC2jcolIt) {
                            size_t commonElSize = elC2jcolIt->first;
                            const map<size_t, size_t>& localIt2keys = elC2jcolIt->second;

                            // the same number of elements is always provided in each call
                            std::vector<std::map<size_t, size_t> > compressPos(commonElSize);
                            std::vector<std::map<size_t, size_t> > resultPos(commonElSize);

                            set<size_t> keys;

                            map<size_t, size_t>::const_iterator lIt2jcolIt;
                            for (lIt2jcolIt = localIt2keys.begin(); lIt2jcolIt != localIt2keys.end(); ++lIt2jcolIt) {
                                size_t localIt = lIt2jcolIt->first;
                                size_t key = lIt2jcolIt->second;

                                keys.insert(key);

                                const std::vector<std::set<size_t> >& origPos = userElLocation.at(key);
                                const set<size_t>& compressed = jcols2e.at(key);

                                set<size_t>::const_iterator itE;
                                size_t e = 0;
                                for (itE = compressed.begin(); itE != compressed.end(); ++itE, e++) {
                                    CPPADCG_ASSERT_UNKNOWN(origPos[*itE].size() == 1);
                                    resultPos[e][localIt] = *origPos[*itE].begin();

                                    compressPos[e][localIt] = *itE;
                                }
                            }

                            ArrayElementGroup* eg = new ArrayElementGroup(keys, commonElSize);
                            data->elCount2elements.m[commonElSize] = eg;

                            for (size_t e = 0; e < commonElSize; e++) {
                                eg->elements[e].resultPattern = IndexPattern::detect(resultPos[e]);
                                eg->elements[e].compressedPattern = IndexPattern::detect(compressPos[e]);
                            }

                        }

                    }

                    // group by number of iterations
                    loopCalls[localit2jcols.size()][loop][group] = data.get();
                    garbage.v.push_back(data.release());
                }
            }
        }

        /**
         * @param loopGroups Used elements from the arrays provided by the group
         *                   function calls (loop->group->{array->{compressed position} })
         * @param nonLoopElements Used elements from non loop function calls
         *                        ([array]{compressed position})
         */
        template<class Base>
        void printForRevUsageFunction(std::ostringstream& out,
                                      const std::string& baseTypeName,
                                      const std::string& modelName,
                                      const std::string& modelFunction,
                                      size_t inLocalSize,
                                      const std::string& localFunction,
                                      const std::string& suffix,
                                      const std::string& keyIndexName,
                                      const IndexDclrOperationNode<Base>& indexIt,
                                      const std::string& resultName,
                                      const std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >& loopGroups,
                                      const std::map<size_t, std::set<size_t> >& nonLoopElements,
                                      const std::map<size_t, std::vector<std::set<size_t> > >& userElLocation,
                                      const std::map<size_t, bool>& jcolOrdered,
                                      void (*generateLocalFunctionName)(std::ostringstream& cache, const std::string& modelName, const LoopModel<Base>& loop, size_t g),
                                      size_t nnz,
                                      size_t maxCompressedSize) {

            using namespace std;

            /**
             * determine to which functions we can provide the Hessian row directly
             * without needing a temporary array (compressed)
             */
            SmartVectorPointer<ArrayGroup> garbage;
            map<size_t, map<LoopModel<Base>*, map<size_t, ArrayGroup*> > > loopCalls;

            /**
             * Determine jrow index patterns and
             * Hessian row start patterns
             */
            determineForRevUsagePatterns(loopGroups, userElLocation, jcolOrdered,
                                         loopCalls, garbage);


            const string& itName = *indexIt.getName();

            string nlRev2Suffix = "noloop_" + suffix;

            CLanguage<Base> langC(baseTypeName);
            string loopFArgs = "inLocal, outLocal, " + langC.getArgumentAtomic();
            string argsDcl = langC.generateDefaultFunctionArgumentsDcl();

            out << "void " << modelFunction << "(" << argsDcl << ") {\n";

            /**
             * Find random index patterns
             */
            set<RandomIndexPattern*> indexRandomPatterns;

            typename map<size_t, map<LoopModel<Base>*, map<size_t, ArrayGroup*> > >::const_iterator itItlg;
            typename map<LoopModel<Base>*, map<size_t, ArrayGroup*> >::const_iterator itlg;
            typename map<size_t, ArrayGroup*>::const_iterator itg;

            for (itItlg = loopCalls.begin(); itItlg != loopCalls.end(); ++itItlg) {

                for (itlg = itItlg->second.begin(); itlg != itItlg->second.end(); ++itlg) {

                    for (itg = itlg->second.begin(); itg != itlg->second.end(); ++itg) {
                        ArrayGroup* group = itg->second;

                        CodeHandler<Base>::findRandomIndexPatterns(group->pattern.get(), indexRandomPatterns);

                        if (group->startLocPattern.get() != NULL) {
                            CodeHandler<Base>::findRandomIndexPatterns(group->startLocPattern.get(), indexRandomPatterns);

                        } else {
                            map<size_t, ArrayElementGroup*>::const_iterator itc;
                            for (itc = group->elCount2elements.m.begin(); itc != group->elCount2elements.m.end(); ++itc) {
                                const ArrayElementGroup* eg = itc->second;

                                for (size_t e = 0; e < eg->elements.size(); e++) {
                                    const ArrayElementCopyPattern& ePos = eg->elements[e];

                                    CodeHandler<Base>::findRandomIndexPatterns(ePos.resultPattern, indexRandomPatterns);
                                    CodeHandler<Base>::findRandomIndexPatterns(ePos.compressedPattern, indexRandomPatterns);
                                }
                            }
                        }
                    }
                }
            }

            /**
             * static variables
             */
            CLanguage<Base>::generateNames4RandomIndexPatterns(indexRandomPatterns);
            CLanguage<Base>::printRandomIndexPatternDeclaration(out, "   ", indexRandomPatterns);

            /**
             * local variables
             */
            out << "   " << baseTypeName << " const * inLocal[" << inLocalSize << "];\n"
                    "   " << baseTypeName << " inLocal1 = 1;\n"
                    "   " << baseTypeName << " * outLocal[1];\n"
                    "   unsigned long " << itName << ";\n"
                    "   unsigned long " << keyIndexName << ";\n"
                    "   unsigned long e;\n";
            CPPADCG_ASSERT_UNKNOWN(itName != "e" && keyIndexName != "e");
            if (maxCompressedSize > 0) {
                out << "   " << baseTypeName << " compressed[" << maxCompressedSize << "];\n";
            }
            out << "   " << baseTypeName << " * " << resultName << " = out[0];\n"
                    "\n"
                    "   inLocal[0] = in[0];\n"
                    "   inLocal[1] = &inLocal1;\n";
            for (size_t j = 2; j < inLocalSize; j++)
                out << "   inLocal[" << j << "] = in[" << (j - 1) << "];\n";

            out << "\n";

            /**
             * zero the output
             */
            out << "   for(e = 0; e < " << nnz << "; e++) " << resultName << "[e] = 0;\n"
                    "\n";

            /**
             * contributions from equations NOT belonging to loops
             * (must come before the loop related values because of the assignments)
             */
            langC.setArgumentIn("inLocal");
            langC.setArgumentOut("outLocal");
            string argsLocal = langC.generateDefaultFunctionArguments();

            bool lastCompressed = false;
            map<size_t, set<size_t> >::const_iterator it;
            for (it = nonLoopElements.begin(); it != nonLoopElements.end(); ++it) {
                size_t index = it->first;
                const set<size_t>& elPos = it->second;
                const std::vector<set<size_t> >& location = userElLocation.at(index);
                CPPADCG_ASSERT_UNKNOWN(elPos.size() <= location.size()); // it can be lower because not all elements have to be assigned
                CPPADCG_ASSERT_UNKNOWN(elPos.size() > 0);
                bool rowOrdered = jcolOrdered.at(index);

                out << "\n";
                if (rowOrdered) {
                    out << "   outLocal[0] = &" << resultName << "[" << *location[0].begin() << "];\n";
                } else if (!lastCompressed) {
                    out << "   outLocal[0] = compressed;\n";
                }
                out << "   " << localFunction << "_" << nlRev2Suffix << index << "(" << argsLocal << ");\n";
                if (!rowOrdered) {
                    for (set<size_t>::const_iterator itEl = elPos.begin(); itEl != elPos.end(); ++itEl) {
                        size_t e = *itEl;
                        out << "   ";
                        set<size_t>::const_iterator itl;
                        for (itl = location[e].begin(); itl != location[e].end(); ++itl) {
                            out << resultName << "[" << (*itl) << "] += compressed[" << e << "];\n";
                        }
                    }
                }
                lastCompressed = !rowOrdered;
            }

            /**
             * loop related values
             */
            for (itItlg = loopCalls.begin(); itItlg != loopCalls.end(); ++itItlg) {
                size_t itCount = itItlg->first;
                if (itCount > 1) {
                    lastCompressed = false;
                    out << "   for(" << itName << " = 0; " << itName << " < " << itCount << "; " << itName << "++) {\n";
                }

                for (itlg = itItlg->second.begin(); itlg != itItlg->second.end(); ++itlg) {
                    LoopModel<Base>& loop = *itlg->first;

                    for (itg = itlg->second.begin(); itg != itlg->second.end(); ++itg) {
                        size_t g = itg->first;
                        ArrayGroup* group = itg->second;

                        const map<size_t, set<size_t> >& key2Compressed = loopGroups.at(&loop).at(g);

                        string indent = itCount == 1 ? "   " : "      "; //indentation

                        if (group->startLocPattern.get() != NULL) {
                            // determine hessRowStart = f(it)
                            out << indent << "outLocal[0] = &" << resultName << "[" << CLanguage<Base>::indexPattern2String(*group->startLocPattern, indexIt) << "];\n";
                        } else {
                            if (!lastCompressed) {
                                out << indent << "outLocal[0] = compressed;\n";
                            }
                            out << indent << "for(e = 0; e < " << maxCompressedSize << "; e++)  compressed[e] = 0;\n";
                        }

                        if (itCount > 1) {
                            out << indent << keyIndexName << " = " << CLanguage<Base>::indexPattern2String(*group->pattern, indexIt) << ";\n";
                            out << indent;
                            (*generateLocalFunctionName)(out, modelName, loop, g);
                            out << "(" << keyIndexName << ", " << loopFArgs << ");\n";
                        } else {
                            size_t key = key2Compressed.begin()->first; // only one jrow
                            out << indent;
                            (*generateLocalFunctionName)(out, modelName, loop, g);
                            out << "(" << key << ", " << loopFArgs << ");\n";
                        }

                        if (group->startLocPattern.get() == NULL) {
                            CPPADCG_ASSERT_UNKNOWN(!group->elCount2elements.m.empty());

                            std::set<size_t> usedIter;

                            // add keys which are never used to usedIter to improve the if/else condition
                            size_t eKey = 0;
                            map<size_t, set<size_t> >::const_iterator itKey = key2Compressed.begin();
                            for (; itKey != key2Compressed.end(); ++itKey) {
                                size_t key = itKey->first;
                                for (size_t k = eKey; k < key; k++) {
                                    usedIter.insert(k);
                                }
                                eKey = key + 1;
                            }

                            bool withIfs = group->elCount2elements.m.size() > 1;
                            map<size_t, ArrayElementGroup*>::const_iterator itc;
                            for (itc = group->elCount2elements.m.begin(); itc != group->elCount2elements.m.end(); ++itc) {
                                const ArrayElementGroup* eg = itc->second;
                                CPPADCG_ASSERT_UNKNOWN(!eg->elements.empty());

                                string indent2 = indent;
                                if (withIfs) {
                                    out << indent;
                                    if (itc != group->elCount2elements.m.begin())
                                        out << "} else ";
                                    if (itc->first != group->elCount2elements.m.rbegin()->first) { // check that it is not the last branch
                                        out << "if(";

                                        size_t maxKey = key2Compressed.rbegin()->first;
                                        std::vector<size_t> info = createIndexConditionExpression(eg->keys, usedIter, maxKey);
                                        CLanguage<Base>::printIndexCondExpr(out, info, keyIndexName);
                                        out << ") ";

                                        usedIter.insert(eg->keys.begin(), eg->keys.end());
                                    }
                                    out << "{\n";
                                    indent2 += "   ";
                                }

                                for (size_t e = 0; e < eg->elements.size(); e++) {
                                    const ArrayElementCopyPattern& ePos = eg->elements[e];

                                    out << indent2 << resultName << "["
                                            << CLanguage<Base>::indexPattern2String(*ePos.resultPattern, indexIt)
                                            << "] += compressed["
                                            << CLanguage<Base>::indexPattern2String(*ePos.compressedPattern, indexIt)
                                            << "];\n";
                                }
                            }

                            if (withIfs) {
                                out << indent << "}\n";
                            }
                        }

                        out << "\n";

                        lastCompressed = group->startLocPattern.get() == NULL;
                    }
                }

                if (itCount > 1) {
                    out << "   }\n";
                }
            }

            out << "\n"
                    "}\n";
        }

        /**
         * 
         * @param elements
         * @param loopGroups Used elements from the arrays provided by the group
         *                   function calls (loop->group->{array->{compressed position} })
         * @param nonLoopElements Used elements from non loop function calls
         *                        ([array]{compressed position})
         * @param functionName
         * @param modelName
         * @param baseTypeName
         * @param suffix
         * @param generateLocalFunctionName
         * @return 
         */
        template<class Base>
        std::string generateGlobalForRevWithLoopsFunctionSource(const std::map<size_t, std::vector<size_t> >& elements,
                                                                const std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >& loopGroups,
                                                                const std::map<size_t, std::set<size_t> >& nonLoopElements,
                                                                const std::string& functionName,
                                                                const std::string& modelName,
                                                                const std::string& baseTypeName,
                                                                const std::string& suffix,
                                                                void (*generateLocalFunctionName)(std::ostringstream& cache, const std::string& modelName, const LoopModel<Base>& loop, size_t g)) {

            using namespace std;

            // functions for each row
            map<size_t, map<LoopModel<Base>*, set<size_t> > > functions;

            typename map<LoopModel<Base>*, map<size_t, map<size_t, set<size_t> > > >::const_iterator itlj1g;
            for (itlj1g = loopGroups.begin(); itlj1g != loopGroups.end(); ++itlj1g) {
                LoopModel<Base>* loop = itlj1g->first;

                map<size_t, map<size_t, set<size_t> > >::const_iterator itg;
                for (itg = itlj1g->second.begin(); itg != itlj1g->second.end(); ++itg) {
                    size_t group = itg->first;
                    const map<size_t, set<size_t> >& jrows = itg->second;

                    for (map<size_t, set<size_t> >::const_iterator itJrow = jrows.begin(); itJrow != jrows.end(); ++itJrow) {
                        functions[itJrow->first][loop].insert(group);
                    }
                }
            }

            /**
             * The function that matches each equation to a directional derivative function
             */
            CLanguage<Base> langC(baseTypeName);
            string argsDcl = langC.generateDefaultFunctionArgumentsDcl();
            string args = langC.generateDefaultFunctionArguments();
            string noLoopFunc = functionName + "_noloop_" + suffix;

            std::ostringstream out;
            out << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n"
                    "\n";
            CLangCompileModelHelper<Base>::generateFunctionDeclarationSource(out, functionName, "noloop_" + suffix, nonLoopElements, argsDcl);
            generateFunctionDeclarationSourceLoopForRev(out, langC, modelName, "j", loopGroups, generateLocalFunctionName);
            out << "\n";
            out << "int " << functionName <<
                    "(unsigned long pos, " << argsDcl << ") {\n"
                    "   \n"
                    "   switch(pos) {\n";
            map<size_t, std::vector<size_t> >::const_iterator it;
            for (it = elements.begin(); it != elements.end(); ++it) {
                size_t jrow = it->first;
                // the size of each sparsity row
                out << "      case " << jrow << ":\n";

                /**
                 * contributions from equations not in loops 
                 * (must come before contributions from loops because of the assignments)
                 */
                map<size_t, set<size_t> >::const_iterator itnl = nonLoopElements.find(jrow);
                if (itnl != nonLoopElements.end()) {
                    out << "         " << noLoopFunc << jrow << "(" << args << ");\n";
                }

                /**
                 * contributions from equations in loops
                 */
                const map<LoopModel<Base>*, set<size_t> >& rowFunctions = functions[jrow];

                typename map<LoopModel<Base>*, set<size_t> >::const_iterator itlg;
                for (itlg = rowFunctions.begin(); itlg != rowFunctions.end(); ++itlg) {
                    LoopModel<Base>* loop = itlg->first;

                    for (set<size_t>::const_iterator itg = itlg->second.begin(); itg != itlg->second.end(); ++itg) {
                        out << "         ";
                        generateLocalFunctionName(out, modelName, *loop, *itg);
                        out << "(" << jrow << ", " << args << ");\n";
                    }
                }

                /**
                 * return all OK
                 */
                out << "         return 0; // done\n";
            }
            out << "      default:\n"
                    "         return 1; // error\n"
                    "   };\n";

            out << "}\n";
            return out.str();
        }

        template<class Base>
        void generateFunctionDeclarationSourceLoopForRev(std::ostringstream& out,
                                                         CLanguage<Base>& langC,
                                                         const std::string& modelName,
                                                         const std::string& keyName,
                                                         const std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >& loopGroups,
                                                         void (*generateLocalFunctionName)(std::ostringstream& cache, const std::string& modelName, const LoopModel<Base>& loop, size_t g)) {

            std::string argsDcl = langC.generateFunctionArgumentsDcl();
            std::string argsDclLoop = "unsigned long " + keyName + ", " + argsDcl;

            typename std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >::const_iterator itlg;
            for (itlg = loopGroups.begin(); itlg != loopGroups.end(); ++itlg) {
                const LoopModel<Base>& loop = *itlg->first;

                typename std::map<size_t, std::map<size_t, std::set<size_t> > >::const_iterator itg;
                for (itg = itlg->second.begin(); itg != itlg->second.end(); ++itg) {
                    size_t group = itg->first;

                    out << "void ";
                    (*generateLocalFunctionName)(out, modelName, loop, group);
                    out << "(" << argsDclLoop << ");\n";
                }
            }
        }

    }
}

#endif
