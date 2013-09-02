#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_REV2_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_REV2_INCLUDED
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

    /***************************************************************************
     *  Utility classes
     **************************************************************************/

    template<class Base>
    class IndexedDependentLoopInfo2 {
    public:
        // {location; origVal}
        typedef std::pair<size_t, CG<Base> > DepLocation;
    public:
        std::map<size_t, DepLocation> iteration2location;
        IndexPattern* pattern;

        inline IndexedDependentLoopInfo2() :
            pattern(NULL) {
        }
    };

    template<class Base>
    class LoopRev2GroupValInfo {
    public:
        CppAD::CG<Base> val;
        size_t location;
        size_t jrow1;
    };

    template<class Base>
    class GroupLoopRev2ColInfo {
    public:
        typedef std::pair<size_t, size_t> TapeVarType;
    public:
        size_t refJrow;
        size_t refIteration;
        std::set<TapeVarType> tapeJ2OrigJ2;
        std::map<size_t, std::map<size_t, std::map<TapeVarType, LoopRev2GroupValInfo<Base> > > > jRow2iterations2Loc;
        std::map<size_t, std::map<size_t, size_t> > jrow2localIt2ModelIt;
    public:

        GroupLoopRev2ColInfo() {
        }

        GroupLoopRev2ColInfo(const std::set<TapeVarType>& vars,
                             const std::map<size_t, std::map<size_t, std::map<TapeVarType, LoopRev2GroupValInfo<Base> > > >& jRowiters) :
            refJrow(jRowiters.begin()->first),
            refIteration(jRowiters.begin()->second.begin()->first),
            tapeJ2OrigJ2(vars),
            jRow2iterations2Loc(jRowiters) {
        }

    };

    /***************************************************************************
     *  Methods related with loop insertion into the operation graph
     **************************************************************************/

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareSparseReverseTwoWithLoops(std::map<std::string, std::string>& sources,
                                                                         const std::map<size_t, std::vector<size_t> >& elements) {
        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        size_t mLoopTotal = 0;
        std::vector<bool> eqLoop(m, false);

        /**
         * loops
         */
        typename std::set<LoopAtomicFun<Base>* >::const_iterator itl;
        for (itl = _loopAtomics.begin(); itl != _loopAtomics.end(); ++itl) {
            LoopAtomicFun<Base>* loop = *itl;

            mLoopTotal += loop->getLoopDependentCount();

            CodeHandler<Base> handler;
            handler.setVerbose(_verbose);

            std::map<size_t, std::vector<std::pair<size_t, CG<Base> > > > hess;

            std::vector<CGBase> tx1v(n, Base(0));

            std::vector<CGBase> tx0(n);
            handler.makeVariables(tx0);
            if (_x.size() > 0) {
                for (size_t i = 0; i < n; i++) {
                    tx0[i].setValue(_x[i]);
                }
            }

            CGBase tx1;
            handler.makeVariable(tx1);
            if (_x.size() > 0) {
                tx1.setValue(Base(1.0));
            }

            /**
             * Generate one function for each dependent variable
             */
            std::map<size_t, std::vector<size_t> >::const_iterator it;
            for (it = elements.begin(); it != elements.end(); ++it) {
                size_t j = it->first;
                const std::vector<size_t>& cols = it->second;
                std::vector<std::pair<size_t, CG<Base> > >& hessRow = hess[j];

                _cache.str("");
                _cache << "model (reverse two, indep " << j << ")";
                const std::string jobName = _cache.str();

                startingGraphCreation(jobName);

                std::vector<CGBase> py(m); // (k+1)*m is not used because we are not interested in all values

                const std::vector<std::vector<LoopPosition> >& eqIndexes = loop->getDependentIndexes();
                for (size_t eq = 0; eq < eqIndexes.size(); eq++) {
                    for (size_t it = 0; it < eqIndexes[eq].size(); it++) {
                        size_t orig = eqIndexes[eq][it].original;
                        eqLoop[orig] = true;

                        handler.makeVariable(py[orig]);
                        if (_x.size() > 0) {
                            py[orig].setValue(Base(1.0));
                        }
                    }
                }

                _fun->Forward(0, tx0);

                tx1v[j] = tx1;
                _fun->Forward(1, tx1v);
                tx1v[j] = Base(0);
                std::vector<CGBase> px = _fun->Reverse(2, py);
                assert(px.size() == 2 * n);

                hessRow.resize(cols.size());
                std::vector<size_t>::const_iterator it2;
                size_t e = 0;
                for (it2 = cols.begin(); it2 != cols.end(); ++it2, e++) {
                    size_t jj = *it2;
                    CGBase& val = px[jj * 2 + 1]; // not interested in all values
                    if (!IdenticalZero(val)) {
                        hessRow[e] = std::make_pair(jj, val);
                    }
                }

                finishedGraphCreation();
            }

            prepareSparseReverseTwoSourcesForLoop(sources, handler, *loop, hess);
        }

        if (mLoopTotal < m) {
            /**
             * equations not in loops
             */
            CodeHandler<Base> handler;
            handler.setVerbose(_verbose);

            std::vector<CGBase> tx1v(n, Base(0));

            /**
             * Generate one function for each dependent variable
             */
            std::map<size_t, std::vector<size_t> >::const_iterator it;
            for (it = elements.begin(); it != elements.end(); ++it) {
                size_t j = it->first;
                const std::vector<size_t>& cols = it->second;

                _cache.str("");
                _cache << "model (reverse two, indep " << j << ")";
                const std::string jobName = _cache.str();

                startingGraphCreation(jobName);

                std::vector<CGBase> tx0(n);
                handler.makeVariables(tx0);
                if (_x.size() > 0) {
                    for (size_t i = 0; i < n; i++) {
                        tx0[i].setValue(_x[i]);
                    }
                }

                CGBase tx1;
                handler.makeVariable(tx1);
                if (_x.size() > 0) {
                    tx1.setValue(Base(1.0));
                }

                std::vector<CGBase> py(m); // (k+1)*m is not used because we are not interested in all values
                for (size_t i = 0; i < m; i++) {
                    if (!eqLoop[i]) {
                        handler.makeVariable(py[i]);
                        if (_x.size() > 0) {
                            py[i].setValue(Base(1.0));
                        }
                    }
                }

                _fun->Forward(0, tx0);

                tx1v[j] = tx1;
                _fun->Forward(1, tx1v);
                tx1v[j] = Base(0);
                std::vector<CGBase> px = _fun->Reverse(2, py);
                assert(px.size() == 2 * n);

                std::vector<CGBase> pxCustom(cols.size());
                std::vector<size_t>::const_iterator it2;
                size_t e = 0;
                size_t nnz = 0;
                for (it2 = cols.begin(); it2 != cols.end(); ++it2, e++) {
                    size_t jj = *it2;
                    CGBase& val = px[jj * 2 + 1]; // not interested in all values
                    if (!IdenticalZero(val)) {
                        pxCustom[e] = val;
                        nnz++;
                    }
                }

                finishedGraphCreation();

                if (nnz > 0) {
                    CLanguage<Base> langC(_baseTypeName);
                    langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
                    _cache.str("");
                    _cache << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO << "_indep" << j;
                    langC.setGenerateFunction(_cache.str());

                    std::ostringstream code;
                    std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("px", "x", "var", "array"));
                    CLangDefaultReverse2VarNameGenerator<Base> nameGenRev2(nameGen.get(), n, 1);

                    handler.generateCode(code, langC, pxCustom, nameGenRev2, _atomicFunctions, jobName);
                }
            }

        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareSparseReverseTwoSourcesForLoop(std::map<std::string, std::string>& sources,
                                                                              CodeHandler<Base>& handler,
                                                                              LoopAtomicFun<Base>& loop,
                                                                              std::map<size_t, std::vector<std::pair<size_t, CG<Base> > > >& rev2) {
        using namespace std;
        using CppAD::vector;

        size_t iterationCount = loop.getIterationCount();
        size_t nnz = _hessSparsity.rows.size();

        if (!loop.isTapeIndependentsFullyDetached2ndOrder()) {
            throw CGException("There are independent variables which appear as indexed and non-indexed for the same equation pattern");
        }

        std::vector<IndexedDependentLoopInfo2<Base> > garbageCollection; ////////// <- rethink this!!!!!!!
        garbageCollection.reserve(nnz);

        // {location; origVal}
        typedef typename IndexedDependentLoopInfo2<Base>::DepLocation DepLocation;


        /**
         * Generate index patterns for the hessian elements resulting from loops
         */
        // loop -> tape independent 1 -> reference orig independent(temporaries only) 1 ->
        //      -> tape independent 2 -> reference orig independent(temporaries only) 2 -> 
        //      -> orig independent 1 -> iteration -> (position, val)
        map<TapeVarType, map<TapeVarType, map<size_t, IndexedDependentLoopInfo2<Base>* > > > hessIndexPatterns;

        map<LoopAtomicFun<Base>*, vector<IndexedDependentLoopInfo2<Base>* > > dependentIndexes;
        // j1 -> iteration -> loop atomic evaluation -> results
        map<TapeVarType, map<size_t, map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> > > > evaluations;

        typename std::map<size_t, std::vector<DepLocation> >::iterator itJRow;
        for (itJRow = rev2.begin(); itJRow != rev2.end(); ++itJRow) {
            size_t j1 = itJRow->first;
            std::vector<DepLocation>& cols = itJRow->second;

            for (size_t e = 0; e < cols.size(); e++) {
                size_t j2 = cols[e].first;
                CGBase& hessVal = cols[e].second;
                if (IdenticalZero(hessVal))
                    continue;

                // loop evaluation -> results
                map<LoopEvaluationOperationNode<Base>*, map<size_t, OperationNode<Base>*> > evals;
                findLoopEvaluations(handler, &loop, hessVal.getOperationNode(), evals); ////////////

                if (evals.empty())
                    continue;

                //size_t iteration = iterationCount;
                if (evals.size() > 1 || evals.begin()->second.size() > 1) {
                    throw CGException("Unable to generate expression for a second order reverse mode element which "
                                      "is either associated with multiple indexed independents "
                                      "or an indepedent with constant index.");
                }

#ifndef NDEBUG
                {
                    LoopEvaluationOperationNode<Base>* loopEvalNode = evals.begin()->first;
                    assert(loopEvalNode->getOperationType() == CGLoopReverseOp && loopEvalNode->getOrder() == 1);
                }
#endif

                /**
                 * If tapeJ1s.size()>1 then the independent is used by several
                 * temporary variables. We know that it is NOT a mix of 
                 * temporaries and indexed independents because this situation
                 * is checked before.
                 */
                set<size_t> tapeJ1s = loop.getIndependentTapeIndexes(j1);
                assert(tapeJ1s.size() >= 1);
                bool isTemporary1 = loop.isTemporary(*tapeJ1s.begin());
                size_t jRef1 = isTemporary1 ? j1 : 0;
                TapeVarType jTape1(*tapeJ1s.begin(), jRef1);

                set<size_t> tapeJ2s = loop.getIndependentTapeIndexes(j2);
                assert(tapeJ2s.size() >= 1); // if >1 then the independent is used by several temporary variables
                bool isTemporary2 = loop.isTemporary(*tapeJ2s.begin());
                size_t jRef2 = isTemporary2 ? j2 : 0;
                TapeVarType jTape2(*tapeJ2s.begin(), jRef2);

                IndexedDependentLoopInfo2<Base>* origRev2El = NULL;

                map<size_t, IndexedDependentLoopInfo2<Base>* >& ref = hessIndexPatterns[jTape1][jTape2];
                if (ref.size() != 0) {
                    typename map<size_t, IndexedDependentLoopInfo2<Base>* >::const_iterator it;
                    it = ref.find(j1);
                    if (it != ref.end()) {
                        origRev2El = it->second;
                    }
                }

                if (origRev2El == NULL) {
                    // the vector will never allocate more space so this is safe:
                    garbageCollection.resize(garbageCollection.size() + 1);
                    ref[j1] = origRev2El = &garbageCollection.back();
                    dependentIndexes[&loop].push_back(origRev2El);
                }

                /**
                 * Determine the iteration
                 */
                size_t iteration;

                bool indexed1 = loop.isIndexedIndependent(*tapeJ1s.begin());
                bool indexed2 = loop.isIndexedIndependent(*tapeJ2s.begin());

                std::set<size_t> iterations;

                if (indexed1 || indexed2) {
                    if (indexed1) {
                        iterations = loop.getIterationsOfIndexedIndep(*tapeJ1s.begin(), j1);
                        assert(!iterations.empty());
                    }

                    if (!indexed1 || (iterations.size() > 1 && indexed2)) {
                        std::set<size_t> iterations2 = loop.getIterationsOfIndexedIndep(*tapeJ2s.begin(), j2);
                        if (iterations.empty()) {
                            iterations = iterations2;
                        } else {
                            std::set<size_t> intersection;
                            std::set_intersection(iterations.begin(), iterations.end(),
                                                  iterations2.begin(), iterations2.end(),
                                                  std::inserter(intersection, intersection.end()));

                            iterations.swap(intersection);
                        }
                    }

                    std::set<size_t>::const_iterator itIt;
                    for (itIt = iterations.begin(); itIt != iterations.end(); ++itIt) {
                        DepLocation& dl = origRev2El->iteration2location[*itIt];
                        dl.first = e; //location
                        dl.second = hessVal; // original value
                    }

                    if (iterations.size() == 1) {
                        iteration = *iterations.begin();
                    } else {
                        if (iterations.find(0) != iterations.end()) {
                            iteration = 0;
                            // must change the argument of the call to the atomic
                            // function so that it only corresponds to a single iteration
                            throw CGException("Not implemented yet");
                        } else {
                            iteration = *iterations.begin();
                        }
                    }

                } else {
                    /**
                     * hessian element relative to two non-indexed independents
                     * (present in all iterations)
                     */
                    iteration = 0; // it is actually present in all iterations!!!
                    for (size_t iter = 0; iter < iterationCount; iter++) {
                        origRev2El->iteration2location[iter] = DepLocation(e, hessVal);
                    }
                    // must change the argument of the call to the atomic
                    // function so that it only corresponds to a single iteration
                    throw CGException("Not implemented yet");
                }

                typename map<LoopEvaluationOperationNode<Base>*, map<size_t, OperationNode<Base>*> >::const_iterator itE;
                for (itE = evals.begin(); itE != evals.end(); ++itE) {
                    LoopEvaluationOperationNode<Base>* loopEvalNode = itE->first;
                    const map<size_t, OperationNode<Base>*>& results = itE->second;

                    size_t result_size = loopEvalNode->getTapeResultSize(); // tape result size
                    vector<OperationNode<Base>*>& allLoopEvalResults = evaluations[jTape1][iteration][loopEvalNode];
                    allLoopEvalResults.resize(result_size);

                    typename map<size_t, OperationNode<Base>*>::const_iterator itR;
                    for (itR = results.begin(); itR != results.end(); ++itR) {
                        size_t tapeIndex = itR->second->getInfo()[1];
                        allLoopEvalResults[tapeIndex] = itR->second;
                    }
                }

            }
        }

        /**
         * Generate index patterns for the dependent variables
         */
        Index indexJrow("jrow");
        Index indexLocalIt("it");
        Index indexLocalItCount("itCount");

        typename map<TapeVarType, map<TapeVarType, map<size_t, IndexedDependentLoopInfo2<Base>* > > >::iterator itJ1;
        for (itJ1 = hessIndexPatterns.begin(); itJ1 != hessIndexPatterns.end(); ++itJ1) {
            const TapeVarType& jTape1 = itJ1->first;

            // iteration -> {(tapeJ2;origJ2)} -> (val, location, jrow1)
            std::map<size_t, map<TapeVarType, LoopRev2GroupValInfo<Base> > > iter2tapeJ2OrigJ2;

            typename map<TapeVarType, map<size_t, IndexedDependentLoopInfo2<Base>* > >::iterator itJ2;
            for (itJ2 = itJ1->second.begin(); itJ2 != itJ1->second.end(); ++itJ2) {
                const TapeVarType& jTape2 = itJ2->first;
                size_t tapeJ2 = jTape2.first;
                size_t origJ2 = jTape2.second;

                typename map<size_t, IndexedDependentLoopInfo2<Base>* >::iterator itRow;
                for (itRow = itJ2->second.begin(); itRow != itJ2->second.end(); ++itRow) {
                    size_t jrow = itRow->first; //j1
                    IndexedDependentLoopInfo2<Base>* idli = itRow->second;

                    std::map<size_t, size_t> iteration2pos;
                    typename std::map<size_t, typename IndexedDependentLoopInfo2<Base>::DepLocation>::const_iterator itloc;
                    for (itloc = idli->iteration2location.begin(); itloc != idli->iteration2location.end(); ++itloc) {
                        size_t iteration = itloc->first;
                        const DepLocation& loc = itloc->second;
                        iteration2pos[iteration] = loc.first; //location
                    }

                    idli->pattern = IndexPattern::detect(LoopAtomicFun<Base>::ITERATION_INDEX, iteration2pos);
                    handler.manageLoopDependentIndexPattern(idli->pattern);

                    /**
                     * 
                     */
                    // also need to save jrow -> position -> iteration
                    for (itloc = idli->iteration2location.begin(); itloc != idli->iteration2location.end(); ++itloc) {
                        size_t iteration = itloc->first;
                        const DepLocation& loc = itloc->second; // (position; value)

                        TapeVarType j2Tape(tapeJ2, origJ2);
                        LoopRev2GroupValInfo<Base>& info = iter2tapeJ2OrigJ2[iteration][j2Tape];
                        info.location = loc.first;
                        info.val = loc.second;
                        info.jrow1 = jrow;
                    }
                }

            }

            /**
             * J2 variable groups
             */
            vector<GroupLoopRev2ColInfo<Base> > groups;

            typename std::map<size_t, std::map<TapeVarType, LoopRev2GroupValInfo<Base> > >::const_iterator itItj2Tape;
            for (itItj2Tape = iter2tapeJ2OrigJ2.begin(); itItj2Tape != iter2tapeJ2OrigJ2.end(); ++itItj2Tape) {
                // find a group with the same j2 variables
                size_t iteration = itItj2Tape->first;
                const std::map<TapeVarType, LoopRev2GroupValInfo<Base> >& tapeJ2OrigJ2 = itItj2Tape->second;
                std::set<TapeVarType> j2Vars;
                typename std::map<TapeVarType, LoopRev2GroupValInfo<Base> >::const_iterator itJ2Tape;
                for (itJ2Tape = tapeJ2OrigJ2.begin(); itJ2Tape != tapeJ2OrigJ2.end(); ++itJ2Tape) {
                    j2Vars.insert(j2Vars.end(), itJ2Tape->first);
                }

                bool found = false;

                for (long g = groups.size() - 1; g >= 0; g--) {
                    GroupLoopRev2ColInfo<Base>& group = groups[g];
                    if (j2Vars.size() == group.tapeJ2OrigJ2.size() &&
                            std::equal(j2Vars.begin(), j2Vars.end(), group.tapeJ2OrigJ2.begin())) {

                        for (itJ2Tape = tapeJ2OrigJ2.begin(); itJ2Tape != tapeJ2OrigJ2.end(); ++itJ2Tape) {
                            // the same combination of variables exists
                            const TapeVarType& j2Tape = itJ2Tape->first;
                            const LoopRev2GroupValInfo<Base>& info = itJ2Tape->second;
                            group.jRow2iterations2Loc[info.jrow1][iteration][j2Tape] = info;
                        }
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    // new combination of variables
                    std::map<size_t, std::map<size_t, std::map<TapeVarType, LoopRev2GroupValInfo<Base> > > > jRow2iterations2Loc;
                    for (itJ2Tape = tapeJ2OrigJ2.begin(); itJ2Tape != tapeJ2OrigJ2.end(); ++itJ2Tape) {
                        const TapeVarType& j2Tape = itJ2Tape->first;
                        const LoopRev2GroupValInfo<Base>& info = itJ2Tape->second;
                        jRow2iterations2Loc[info.jrow1][iteration][j2Tape] = info;
                    }

                    groups.push_back(GroupLoopRev2ColInfo<Base>(j2Vars, jRow2iterations2Loc));
                }

            }

            /**
             * Generate index patterns
             */
            long g_size = groups.size();
            for (long g = 0; g < g_size; g++) {
                GroupLoopRev2ColInfo<Base>& group = groups[g];

                /**
                 * Model index pattern
                 * 
                 * detect the index pattern for the model iterations
                 * based on jrow and the local loop iteration
                 */
                typename std::map<size_t, std::map<size_t, std::map<TapeVarType, LoopRev2GroupValInfo<Base> > > >::const_iterator itJrow;
                for (itJrow = group.jRow2iterations2Loc.begin(); itJrow != group.jRow2iterations2Loc.end(); ++itJrow) {
                    size_t jrow = itJrow->first;

                    std::map<size_t, size_t>& localIt2ModelIt = group.jrow2localIt2ModelIt[jrow];
                    typename std::map<size_t, std::map<TapeVarType, LoopRev2GroupValInfo<Base> > >::const_iterator itIt;
                    size_t localIt = 0;
                    for (itIt = itJrow->second.begin(); itIt != itJrow->second.end(); ++itIt, localIt++) {
                        size_t modelIt = itIt->first;
                        localIt2ModelIt[localIt] = modelIt;
                    }
                }

                /**
                 * try to fit a combination of two patterns:
                 *  j = fStart(jrow) + flit(lit);
                 */
                std::auto_ptr<IndexPattern> itPattern(Plane2DIndexPattern::detectPlane2D(indexJrow, indexLocalIt, group.jrow2localIt2ModelIt));

                if (itPattern.get() == NULL) {
                    // did not match!
                    throw CGException("Not implemented yet");
                }

                /**
                 * Local iteration count pattern
                 */
                std::auto_ptr<IndexPattern> indexLocalItCountPattern;
                if (group.jrow2localIt2ModelIt.size() > 0) {
                    std::map<size_t, size_t> jrow2litCount;

                    std::map<size_t, std::map<size_t, size_t> >::const_iterator itJrow2Count;

                    for (itJrow2Count = group.jrow2localIt2ModelIt.begin(); itJrow2Count != group.jrow2localIt2ModelIt.end(); ++itJrow2Count) {
                        size_t jrow = itJrow2Count->first;
                        jrow2litCount[jrow] = itJrow2Count->second.size();
                    }

                    indexLocalItCountPattern.reset(IndexPattern::detect(indexLocalItCount, jrow2litCount));
                }

                /**
                 * Compressed row element position index pattern
                 * 
                 */
                std::map<TapeVarType, std::map<size_t, std::map<size_t, size_t> > > el2jrow2localIt2Pos;

                for (itJrow = group.jRow2iterations2Loc.begin(); itJrow != group.jRow2iterations2Loc.end(); ++itJrow) {
                    size_t jrow = itJrow->first;

                    typename std::map<size_t, std::map<TapeVarType, LoopRev2GroupValInfo<Base> > >::const_iterator itIt;
                    size_t localIt = 0;
                    for (itIt = itJrow->second.begin(); itIt != itJrow->second.end(); ++itIt, localIt++) {

                        typename std::map<TapeVarType, LoopRev2GroupValInfo<Base> >::const_iterator itVar;
                        for (itVar = itIt->second.begin(); itVar != itIt->second.end(); ++itVar) {
                            const TapeVarType& el = itVar->first;
                            el2jrow2localIt2Pos[el][jrow][localIt] = itVar->second.location;
                        }

                    }
                }

                SmartMapValuePointer<TapeVarType, Plane2DIndexPattern> loopDepIndexes;
                std::map<TapeVarType, std::map<size_t, std::map<size_t, size_t> > >::const_iterator itEl;
                for (itEl = el2jrow2localIt2Pos.begin(); itEl != el2jrow2localIt2Pos.end(); ++itEl) {
                    const TapeVarType& el = itEl->first;
                    const std::map<size_t, std::map<size_t, size_t> >& jrow2localIt2Pos = itEl->second;

                    Plane2DIndexPattern * posPattern2D(Plane2DIndexPattern::detectPlane2D(indexJrow, indexLocalIt, jrow2localIt2Pos));

                    if (posPattern2D == NULL) {
                        // did not match!
                        throw CGException("Not implemented yet");
                    }
                    loopDepIndexes.m[el] = posPattern2D;
                }

                _cache.str("");
                _cache << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO <<
                        "_loop" << loop.getLoopId() <<
                        "_indep" << jTape1.first << "_" << jTape1.second <<
                        "_g" << g;
                std::string functionName = _cache.str();

                _cache.str("");
                _cache << "model (reverse two, loop " << loop.getLoopId() << ", indep " << jTape1.first << "_" << jTape1.second << ", group " << g << ")";
                std::string jobName = _cache.str();

                std::string source = generateSparseReverseTwoWithLoopsVarGroupSource(functionName,
                                                                                     jobName,
                                                                                     loop, handler,
                                                                                     jTape1, group,
                                                                                     indexJrow,
                                                                                     indexLocalIt,
                                                                                     indexLocalItCount,
                                                                                     *itPattern.get(),
                                                                                     indexLocalItCountPattern.get(),
                                                                                     loopDepIndexes.m,
                                                                                     evaluations[jTape1][group.refIteration]);

                sources[functionName + ".c"] = source;
            }
        }

        //size_t assignOrAdd = 1; // add
        //prepareLoops(handler, hess, evaluations, dependentIndexes, assignOrAdd);
    }

    template<class Base>
    std::string CLangCompileModelHelper<Base>::generateSparseReverseTwoWithLoopsVarGroupSource(const std::string& functionName,
                                                                                               const std::string& jobName,
                                                                                               LoopAtomicFun<Base>& loop,
                                                                                               CodeHandler<Base>& handler,
                                                                                               const TapeVarType& jTape1,
                                                                                               const GroupLoopRev2ColInfo<Base>& group,
                                                                                               const Index& indexJrow,
                                                                                               const Index& indexLocalIt,
                                                                                               const Index& indexLocalItCount,
                                                                                               const IndexPattern& itPattern,
                                                                                               const IndexPattern* itCountPattern,
                                                                                               const std::map<TapeVarType, Plane2DIndexPattern*>& loopDepIndexes,
                                                                                               std::map<LoopEvaluationOperationNode<Base>*, vector<OperationNode<Base>*> >& evaluations) {
        size_t n = _fun->Domain();

        /**
         * determine if a loop should be created
         */
        OperationNode<Base>* loopStart = NULL;

        std::map<size_t, std::set<size_t> > localIterCount2Jrows;

        std::map<size_t, std::map<size_t, size_t> >::const_iterator itJrow2localit;
        for (itJrow2localit = group.jrow2localIt2ModelIt.begin(); itJrow2localit != group.jrow2localIt2ModelIt.end(); ++itJrow2localit) {
            size_t jrow = itJrow2localit->first;
            size_t itCount = itJrow2localit->second.size();
            localIterCount2Jrows[itCount].insert(jrow);
        }

        bool createsLoop = localIterCount2Jrows.size() != 1 || // is there a different number of iterations?
                localIterCount2Jrows.begin()->first != 1; // is there always just on iteration?

        IndexOperationNode<Base> jrowIndexOp(indexJrow);

        std::auto_ptr<LoopNodeInfo<Base> > loopNodeInfo;
        std::auto_ptr<IndexOperationNode<Base> > localIterIndexOp;
        std::auto_ptr<IndexOperationNode<Base> > localIterCountIndexOp;
        std::auto_ptr<IndexAssignOperationNode<Base> > itCountAssignOp;
        if (createsLoop) {
            if (itCountPattern != NULL) {
                itCountAssignOp.reset(new IndexAssignOperationNode<Base>(indexLocalItCount, *itCountPattern, jrowIndexOp));
                localIterCountIndexOp.reset(new IndexOperationNode<Base>(indexLocalItCount, *itCountAssignOp.get()));
                loopNodeInfo.reset(new CustomLoopNodeInfo<Base>(indexLocalIt, *localIterCountIndexOp.get()));
            } else {
                size_t itCount = group.jrow2localIt2ModelIt.begin()->second.size();
                loopNodeInfo.reset(new CustomLoopNodeInfo<Base>(indexLocalIt, itCount));
            }

            loopStart = new LoopStartOperationNode<Base>(*loopNodeInfo.get());
            handler.manageOperationNodeMemory(loopStart);

            localIterIndexOp.reset(new IndexOperationNode<Base>(indexLocalIt, *loopStart));
        }


        IndexAssignOperationNode<Base> iterationIndexPatternOp(LoopAtomicFun<Base>::ITERATION_INDEX, itPattern, &jrowIndexOp, localIterIndexOp.get());
        IndexOperationNode<Base> iterationIndexOp(LoopAtomicFun<Base>::ITERATION_INDEX, iterationIndexPatternOp);

        /**
         * generate the expressions inside the loop
         */
        vector<OperationNode<Base>* > indexedIndependents;
        vector<OperationNode<Base>* > indexedIndependents2;

        replaceAtomicLoopWithExpression(handler, loop, iterationIndexOp, evaluations, indexedIndependents, indexedIndependents2);

        /**
         * 
         */
        const std::map<TapeVarType, LoopRev2GroupValInfo<Base> >& refVals = group.jRow2iterations2Loc.at(group.refJrow).at(group.refIteration);
        std::vector<CGBase> pxCustom(refVals.size());

        size_t assignOrAdd = 1; //add

        std::map<TapeVarType, Plane2DIndexPattern*>::const_iterator itDepIndex;
        /**
         * make the loop end
         */
        if (createsLoop) {
            vector<std::pair<CGBase, IndexPattern*> > indexedLoopResults(refVals.size());

            size_t i = 0;
            for (itDepIndex = loopDepIndexes.begin(); itDepIndex != loopDepIndexes.end(); ++itDepIndex, i++) {
                const CGBase& val = refVals.at(itDepIndex->first).val;
                indexedLoopResults[i] = std::make_pair(val, itDepIndex->second);
            }

            OperationNode<Base>* loopEnd = createLoopEnd(handler, indexedLoopResults, *loopNodeInfo.get(), assignOrAdd);

            /**
             * change the dependents (must depend directly on the loop)
             */
            for (size_t i = 0; i < pxCustom.size(); i++) {
                pxCustom[i] = handler.createCG(Argument<Base>(*loopEnd));
            }

            /**
             * move no-nindexed expressions outside loop
             */
            moveNonIndexedOutsideLoop(*loopStart, *loopEnd);
        } else {
            /**
             * No loop
             */
            std::vector<Argument<Base> > indexedArgs(1);
            std::vector<size_t> info(2);

            size_t i = 0;
            for (itDepIndex = loopDepIndexes.begin(); itDepIndex != loopDepIndexes.end(); ++itDepIndex, i++) {
                const CGBase& val = refVals.at(itDepIndex->first).val;

                indexedArgs[0] = asArgument(val); // indexed expression
                info[0] = handler.addLoopDependentIndexPattern(*itDepIndex->second); // dependent index pattern location
                info[1] = assignOrAdd;

                OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, info, indexedArgs);
                handler.manageOperationNodeMemory(yIndexed);

                pxCustom[i] = handler.createCG(Argument<Base>(*yIndexed));
            }
        }

        CLanguage<Base> langC(_baseTypeName);
        langC.setCurrentLoop(&loop);

        _cache.str("");
        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("px", "x", "var", "array"));
        CLangDefaultReverse2VarNameGenerator<Base> nameGenRev2(nameGen.get(), n, 1);

        /**
         * Generate the source code inside the loop
         */
        handler.generateCode(code, langC, pxCustom, nameGenRev2, _atomicFunctions, jobName);

        /**
         * 
         */
        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                "#include <math.h>\n"
                "\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n"
                "\n"
                "void " << functionName << "(unsigned long jrow, " << argsDcl << ") {\n";
        nameGenRev2.customFunctionVariableDeclarations(_cache);
        _cache << langC.generateIndependentVariableDeclaration() << "\n";
        _cache << langC.generateDependentVariableDeclaration() << "\n";
        _cache << langC.generateTemporaryVariableDeclaration(true) << "\n";
        nameGenRev2.prepareCustomFunctionVariables(_cache);

        // code inside the loop
        _cache << code.str();

        if (createsLoop) {
            _cache << "}\n";
        }

        nameGenRev2.finalizeCustomFunctionVariables(_cache);
        _cache << "}\n\n";

        std::string source = _cache.str();
        _cache.str("");
        return source;
    }

    /**
     * 
     * @param modelFunction
     * @param functionRev2
     * @param rev2Suffix
     * @param userHessElLocation maps each element to its position in the user hessian
     * @param elements  elements[var]{var}
     * @return 
     */
    template<class Base>
    std::string CLangCompileModelHelper<Base>::generateSparseHessianSourceFromRev2WithLoops(const std::string& modelFunction,
                                                                                            const std::string& functionRev2,
                                                                                            const std::string& rev2Suffix,
                                                                                            const std::map<size_t, std::vector<std::set<size_t> > >& userHessElLocation,
                                                                                            const std::map<size_t, std::vector<size_t> >& elements,
                                                                                            size_t maxCompressedSize) {

        //size_t m = _fun->Range();
        //size_t n = _fun->Domain();
        using namespace std;
        /*
                CLanguage<Base> langC(_baseTypeName);
                std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();

                _cache.str("");
                _cache << "#include <stdlib.h>\n"
                        << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n\n";
                generateFunctionDeclarationSource(_cache, functionRev2, rev2Suffix, elements, argsDcl);
                _cache << "\n"
                        "void " << modelFunction << "(" << argsDcl << ") {\n"
                        "   " << _baseTypeName << " const * inLocal[3];\n"
                        "   " << _baseTypeName << " inLocal1 = 1;\n"
                        "   " << _baseTypeName << " * outLocal[1];\n"
                        "   " << _baseTypeName << " compressed[" << maxCompressedSize << "];\n"
                        "   " << _baseTypeName << " * hess = out[0];\n"
                        "\n"
                        "   inLocal[0] = in[0];\n"
                        "   inLocal[1] = &inLocal1;\n"
                        "   inLocal[2] = in[1];\n"
                        "   outLocal[0] = compressed;";

                langC.setArgumentIn("inLocal");
                langC.setArgumentOut("outLocal");
                std::string argsLocal = langC.generateDefaultFunctionArguments();

                for (loops) {
                    for (it = elements.begin(); it != elements.end(); ++it) {
                        size_t index = it->first; ////// check if there are elements for the loop
                        const std::vector<size_t>& els = it->second;
                        const std::vector<set<size_t> >& location = userHessElLocation.at(index);
                        assert(els.size() == location.size());





                        _cache << "\n"
                                "   " << functionRev2 << "_" << rev2Suffix << index << "(" << argsLocal << ");\n";
                        for (size_t e = 0; e < els.size(); e++) {
                            _cache << "   ";
                            set<size_t>::const_iterator itl;
                            for (itl = location[e].begin(); itl != location[e].end(); ++itl) {
                                _cache << "hess[" << (*itl) << "] = ";
                            }
                            _cache << "compressed[" << e << "];\n";
                        }
                    }
                }


                 // hessian elements comming from equations not in loops

                for (it = elements.begin(); it != elements.end(); ++it) {
                    size_t index = it->first; ////// check if there are elements for the loop
                    const std::vector<size_t>& els = it->second;
                    const std::vector<set<size_t> >& location = userHessElLocation.at(index);
                    assert(els.size() == location.size());

                    _cache << "\n"
                            "   " << functionRev2 << "_" << rev2Suffix << index << "(" << argsLocal << ");\n";
                    for (size_t e = 0; e < els.size(); e++) {
                        _cache << "   ";
                        set<size_t>::const_iterator itl;
                        for (itl = location[e].begin(); itl != location[e].end(); ++itl) {

                            _cache << "hess[" << (*itl) << "] += "; // notice the +=
                        }
                        _cache << "compressed[" << e << "];\n";
                    }
                }

                _cache << "\n"
                        "}\n";


                _cache.str("");
                std::string source = _cache.str();
                _cache.str("");
                return source;
         */
        return ""; //TODO
    }

}

#endif
