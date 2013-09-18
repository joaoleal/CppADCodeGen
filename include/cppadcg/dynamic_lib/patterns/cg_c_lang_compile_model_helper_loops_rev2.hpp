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
    class LoopRev2GroupValInfo {
    public:
        size_t compressedLoc;
        size_t jrow1;
        CppAD::CG<Base> val;
    };

    template<class Base>
    class LoopRev2ValInfo {
    public:

        size_t compressedLoc;
        size_t jcol2;
        CppAD::CG<Base> val;

        LoopRev2ValInfo() {
        }

        LoopRev2ValInfo(size_t e, size_t orig2, const CppAD::CG<Base>& v) :
            compressedLoc(e),
            jcol2(orig2),
            val(v) {
        }
    };

    template<class Base>
    class IndexedDependentLoopInfo2 {
    public:
        std::map<size_t, LoopRev2ValInfo<Base> > iteration2location;
        IndexPattern* pattern;

        inline IndexedDependentLoopInfo2() :
            pattern(NULL) {
        }
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
        IndexPattern* jrowPattern;
        IndexPattern* hessStartLocPattern;
    public:

        GroupLoopRev2ColInfo() {
        }

        GroupLoopRev2ColInfo(const std::set<TapeVarType>& vars,
                             const std::map<size_t, std::map<size_t, std::map<TapeVarType, LoopRev2GroupValInfo<Base> > > >& jRowiters) :
            refJrow(jRowiters.begin()->first),
            refIteration(jRowiters.begin()->second.begin()->first),
            tapeJ2OrigJ2(vars),
            jRow2iterations2Loc(jRowiters),
            jrowPattern(NULL),
            hessStartLocPattern(NULL) {
        }

        virtual ~GroupLoopRev2ColInfo() {
            delete jrowPattern;
            delete hessStartLocPattern;
        }
    };

    /***************************************************************************
     *  Methods related with loop insertion into the operation graph
     **************************************************************************/

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareSparseReverseTwoWithLoops(std::map<std::string, std::string>& sources,
                                                                         const std::map<size_t, std::vector<size_t> >& elements) {
        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        size_t mLoopTotal = 0;
        std::vector<bool> eqLoop(m, false);

        typedef std::pair<size_t, size_t> Compressed2JColType;

        /***********************************************************************
         * loops
         **********************************************************************/
        typename std::set<LoopModel<Base>* >::const_iterator itl;
        for (itl = _loopTapes.begin(); itl != _loopTapes.end(); ++itl) {
            LoopModel<Base>* loop = *itl;

            const std::vector<std::vector<LoopPosition> >& eqIndexes = loop->getDependentIndexes();

            //mLoopTotal += loop->getLoopDependentCount();

            SparsitySetType loopHessSparsity(n);

            for (size_t eq = 0; eq < eqIndexes.size(); eq++) {
                for (size_t it = 0; it < eqIndexes[eq].size(); it++) {
                    size_t orig = eqIndexes[eq][it].original;
                    eqLoop[orig] = true;
                    addMatrixSparsity(_hessSparsities[orig].sparsity, loopHessSparsity);
                }
            }

            /**
             * determine the sparsity using only the loop equations
             */
            // elements[var]{vars}
            std::map<size_t, std::vector<Compressed2JColType> > loopElements;

            std::map<size_t, std::vector<size_t> >::const_iterator itJrow2jcols;
            for (itJrow2jcols = elements.begin(); itJrow2jcols != elements.end(); ++itJrow2jcols) {
                size_t jrow = itJrow2jcols->first;
                const std::vector<size_t>& jcols = itJrow2jcols->second;
                const std::set<size_t>& row = loopHessSparsity[jrow];

                for (size_t e = 0; e < jcols.size(); e++) {
                    if (row.find(jcols[e]) != row.end()) {
                        loopElements[jrow].push_back(Compressed2JColType(e, jcols[e]));
                    }
                }
            }

            if (loopElements.empty())
                continue; // nothing to do

            CodeHandler<Base> handler;
            handler.setVerbose(_verbose);

            std::map<size_t, std::vector<LoopRev2ValInfo<Base> > > hess;

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
             * 
             */
            std::map<size_t, std::vector<Compressed2JColType> >::const_iterator it;
            for (it = loopElements.begin(); it != loopElements.end(); ++it) {
                size_t j = it->first;
                const std::vector<Compressed2JColType>& cols = it->second;
                std::vector<LoopRev2ValInfo<Base> >& hessRow = hess[j];

                _cache.str("");
                _cache << "model (reverse two, indep " << j << ")";
                const std::string jobName = _cache.str();

                startingGraphCreation(jobName);

                std::vector<CGBase> py(m); // (k+1)*m is not used because we are not interested in all values

                for (size_t eq = 0; eq < eqIndexes.size(); eq++) {
                    for (size_t it = 0; it < eqIndexes[eq].size(); it++) {
                        size_t orig = eqIndexes[eq][it].original;

                        handler.makeVariable(py[orig]);
                        if (_x.size() > 0) {
                            py[orig].setValue(Base(1.0));
                        }
                    }
                }

                _fun.Forward(0, tx0);

                tx1v[j] = tx1;
                _fun.Forward(1, tx1v);
                tx1v[j] = Base(0);
                std::vector<CGBase> px = _fun.Reverse(2, py);
                assert(px.size() == 2 * n);

                hessRow.resize(cols.size());

                std::vector<Compressed2JColType>::const_iterator it2;
                for (it2 = cols.begin(); it2 != cols.end(); ++it2) {
                    size_t e = it2->first;
                    size_t jj = it2->second;

                    CGBase& val = px[jj * 2 + 1]; // not interested in all values
                    if (!IdenticalZero(val)) {
                        hessRow.push_back(LoopRev2ValInfo<Base>(e, jj, val));
                    }
                }

                finishedGraphCreation();
            }

            prepareSparseReverseTwoSourcesForLoop(sources, handler, *loop, hess, tx1);
        }

        if (mLoopTotal < m) {
            /*******************************************************************
             * equations NOT in loops
             ******************************************************************/

            SparsitySetType nonLoopHessSparsity(n);

            for (size_t i = 0; i < m; i++) {
                if (!eqLoop[i]) {
                    addMatrixSparsity(_hessSparsities[i].sparsity, nonLoopHessSparsity);
                }
            }

            /**
             * determine the sparsity using only the equations not in loops
             */
            std::map<size_t, std::vector<size_t> >::const_iterator itJrow2jcols;
            for (itJrow2jcols = elements.begin(); itJrow2jcols != elements.end(); ++itJrow2jcols) {
                size_t jrow = itJrow2jcols->first;
                const std::vector<size_t>& jcols = itJrow2jcols->second;
                const std::set<size_t>& row = nonLoopHessSparsity[jrow];

                for (size_t e = 0; e < jcols.size(); e++) {
                    if (row.find(jcols[e]) != row.end()) {
                        _nonLoopRev2Elements[jrow].push_back(Compressed2JColType(e, jcols[e]));
                    }
                }
            }

            std::vector<CGBase> tx1v(n, Base(0));

            /**
             * Generate one function for each dependent variable
             */
            std::map<size_t, std::vector<Compressed2JColType> >::const_iterator it;
            for (it = _nonLoopRev2Elements.begin(); it != _nonLoopRev2Elements.end(); ++it) {
                size_t j = it->first;
                const std::vector<Compressed2JColType>& cols = it->second;

                _cache.str("");
                _cache << "model (reverse two, indep " << j << ")";
                const std::string jobName = _cache.str();

                startingGraphCreation(jobName);

                CodeHandler<Base> handler;
                handler.setVerbose(_verbose);

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
                handler.makeVariables(py);

                std::vector<CGBase> pyNoLoop(m);
                for (size_t i = 0; i < m; i++) {
                    if (!eqLoop[i]) {
                        pyNoLoop[i] = py[i];
                        if (_x.size() > 0) {
                            pyNoLoop[i].setValue(Base(1.0));
                        }
                    }
                }

                _fun.Forward(0, tx0);

                tx1v[j] = tx1;
                _fun.Forward(1, tx1v);
                tx1v[j] = Base(0);
                std::vector<CGBase> px = _fun.Reverse(2, pyNoLoop);
                assert(px.size() == 2 * n);

                vector<CGBase> pxCustom(cols.size());
                std::vector<Compressed2JColType>::const_iterator it2;
                for (it2 = cols.begin(); it2 != cols.end(); ++it2) {
                    size_t e = it2->first;
                    size_t jj = it2->second;
                    CGBase& val = px[jj * 2 + 1]; // not interested in all values
                    if (!IdenticalZero(val)) {
                        if (pxCustom.size() <= e) {
                            pxCustom.resize(e + 1);
                        }
                        pxCustom[e] = val;
                    }
                }

                finishedGraphCreation();

                CLanguage<Base> langC(_baseTypeName);
                langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
                _cache.str("");
                _cache << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO << "_noloop_indep" << j;
                langC.setGenerateFunction(_cache.str());

                std::ostringstream code;
                std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("px", "x", "var", "array"));
                CLangDefaultReverse2VarNameGenerator<Base> nameGenRev2(nameGen.get(), n, 1);

                handler.generateCode(code, langC, pxCustom, nameGenRev2, _atomicFunctions, jobName);
            }

        }

        generateGlobalReverseTwoWithLoopsFunctionSource(elements, sources);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::prepareSparseReverseTwoSourcesForLoop(std::map<std::string, std::string>& sources,
                                                                              CodeHandler<Base>& handler,
                                                                              LoopModel<Base>& loop,
                                                                              std::map<size_t, std::vector<LoopRev2ValInfo<Base> > >& rev2,
                                                                              const CGBase& tx1) {
        using namespace std;
        using CppAD::vector;

        //size_t iterationCount = loop.getIterationCount();
        size_t nnz = _hessSparsity.rows.size();

        std::vector<IndexedDependentLoopInfo2<Base> > garbageCollection; ////////// <- rethink this!!!!!!!
        garbageCollection.reserve(nnz);

        /**
         * Generate index patterns for the hessian elements resulting from loops
         */
        // loop -> tape independent 1 -> reference orig independent(temporaries only) 1 ->
        //      -> tape independent 2 -> reference orig independent(temporaries only) 2 -> 
        //      -> orig independent 1 -> iteration -> (position, val)
        map<TapeVarType, map<TapeVarType, map<size_t, IndexedDependentLoopInfo2<Base>* > > > hessIndexPatterns;

        map<LoopModel<Base>*, vector<IndexedDependentLoopInfo2<Base>* > > dependentIndexes;

        typename std::map<size_t, std::vector<LoopRev2ValInfo<Base> > >::iterator itJRow;


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
                    typename std::map<size_t, LoopRev2ValInfo<Base> >::const_iterator itloc;
                    for (itloc = idli->iteration2location.begin(); itloc != idli->iteration2location.end(); ++itloc) {
                        size_t iteration = itloc->first;
                        const LoopRev2ValInfo<Base>& loc = itloc->second;
                        iteration2pos[iteration] = loc.compressedLoc; //location
                    }

                    idli->pattern = IndexPattern::detect(LoopModel<Base>::ITERATION_INDEX, iteration2pos);
                    handler.manageLoopDependentIndexPattern(idli->pattern);

                    /**
                     * 
                     */
                    // also need to save jrow -> position -> iteration
                    for (itloc = idli->iteration2location.begin(); itloc != idli->iteration2location.end(); ++itloc) {
                        size_t iteration = itloc->first;
                        const LoopRev2ValInfo<Base>& loc = itloc->second; // (position; value)

                        TapeVarType j2Tape(tapeJ2, origJ2);
                        LoopRev2GroupValInfo<Base>& info = iter2tapeJ2OrigJ2[iteration][j2Tape];
                        info.compressedLoc = loc.compressedLoc;
                        info.val = loc.val;
                        info.jrow1 = jrow;
                    }
                }

            }

            /**
             * J2 variable groups
             */
            vector<GroupLoopRev2ColInfo<Base>* >& groups = _loopRev2Groups[&loop][jTape1];

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
                    GroupLoopRev2ColInfo<Base>& group = *groups[g];
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

                    groups.push_back(new GroupLoopRev2ColInfo<Base>(j2Vars, jRow2iterations2Loc));
                }

            }

            /**
             * Generate index patterns
             */
            long g_size = groups.size();
            for (long g = 0; g < g_size; g++) {
                GroupLoopRev2ColInfo<Base>& group = *groups[g];

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

                    indexLocalItCountPattern.reset(IndexPattern::detect(indexJrow, jrow2litCount));
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
                        //size_t iteration = itIt->first; // model iteration

                        typename std::map<TapeVarType, LoopRev2GroupValInfo<Base> >::const_iterator itVar;
                        for (itVar = itIt->second.begin(); itVar != itIt->second.end(); ++itVar) {
                            const TapeVarType& el = itVar->first;
                            el2jrow2localIt2Pos[el][jrow][localIt] = itVar->second.compressedLoc;
                        }

                    }
                }

                SmartMapValuePointer<TapeVarType, Plane2DIndexPattern> loopDepIndexes;
                std::map<TapeVarType, std::map<size_t, std::map<size_t, size_t> > >::const_iterator itEl;
                for (itEl = el2jrow2localIt2Pos.begin(); itEl != el2jrow2localIt2Pos.end(); ++itEl) {
                    const TapeVarType& el = itEl->first;
                    const std::map<size_t, std::map<size_t, size_t> >& jrow2localIt2Pos = itEl->second;

                    Plane2DIndexPattern* posPattern2D = Plane2DIndexPattern::detectPlane2D(indexJrow, indexLocalIt, jrow2localIt2Pos);

                    if (posPattern2D == NULL) {
                        // did not match!
                        throw CGException("Not implemented yet");
                    }
                    loopDepIndexes.m[el] = posPattern2D;
                }

                _cache.str("");
                generateFunctionNameLoopRev2(_cache, loop, jTape1, g);
                std::string functionName = _cache.str();

                _cache.str("");
                _cache << "model (reverse two, loop " << loop.getLoopId() << ", indep " << jTape1.first << "_" << jTape1.second << ", group " << g << ")";
                std::string jobName = _cache.str();
                /*
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
                                                                                                     evaluations[jTape1][group.refIteration],
                                                                                                     tx1);
                 
                sources[functionName + ".c"] = source;
                 * */
            }
        }

    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateGlobalReverseTwoWithLoopsFunctionSource(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                        std::map<std::string, std::string>& sources) {

        // functions for each row
        std::map<size_t, std::map<LoopModel<Base>*, std::map<TapeVarType, std::map<size_t, GroupLoopRev2ColInfo<Base>* > > > > functions;

        typename std::map<LoopModel<Base>*, std::map<TapeVarType, vector<GroupLoopRev2ColInfo<Base>* > > >::const_iterator itlj1g;
        for (itlj1g = _loopRev2Groups.begin(); itlj1g != _loopRev2Groups.end(); ++itlj1g) {
            LoopModel<Base>* loop = itlj1g->first;

            typename std::map<TapeVarType, vector<GroupLoopRev2ColInfo<Base>* > >::const_iterator itj1g;
            for (itj1g = itlj1g->second.begin(); itj1g != itlj1g->second.end(); ++itj1g) {
                const TapeVarType& jTape1 = itj1g->first;

                const vector<GroupLoopRev2ColInfo<Base>* >& groups = itj1g->second;
                for (size_t g = 0; g < groups.size(); g++) {
                    GroupLoopRev2ColInfo<Base>* group = groups[g];

                    std::map<size_t, std::map<size_t, size_t> >::const_iterator itJrow;
                    for (itJrow = group->jrow2localIt2ModelIt.begin(); itJrow != group->jrow2localIt2ModelIt.end(); ++itJrow) {
                        size_t jrow = itJrow->first;
                        functions[jrow][loop][jTape1][g] = group;
                    }

                }
            }
        }

        /**
         * The function that matches each equation to a directional derivative function
         */
        CLanguage<Base> langC(_baseTypeName);
        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();
        std::string args = langC.generateDefaultFunctionArguments();
        std::string functionRev2 = _name + "_" + FUNCTION_SPARSE_REVERSE_TWO;
        std::string noLoopFunc = functionRev2 + "_noloop_indep";

        _cache.str("");
        _cache << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n"
                "\n";
        generateFunctionDeclarationSource(_cache, functionRev2, "noloop_indep", _nonLoopRev2Elements, argsDcl);
        generateFunctionDeclarationSourceLoopRev2(_cache, langC);
        _cache << "\n";
        _cache << "int " << functionRev2 <<
                "(unsigned long pos, " << argsDcl << ") {\n"
                "   \n"
                "   switch(pos) {\n";
        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t jrow = it->first;
            // the size of each sparsity row
            _cache << "      case " << jrow << ":\n";

            /**
             * contributions from equations not in loops 
             * (must come before contributions from loops because of the assigments)
             */
            std::map<size_t, std::vector<Compressed2JColType> >::const_iterator itnl = _nonLoopRev2Elements.find(jrow);
            if (itnl != _nonLoopRev2Elements.end()) {
                _cache << "         " << noLoopFunc << jrow << "(" << args << ");\n";
            }

            /**
             * contributions from equations in loops
             */
            const std::map<LoopModel<Base>*, std::map<TapeVarType, std::map<size_t, GroupLoopRev2ColInfo<Base>* > > >& rowFunctions = functions[jrow];

            typename std::map<LoopModel<Base>*, std::map<TapeVarType, std::map<size_t, GroupLoopRev2ColInfo<Base>* > > >::const_iterator itlj1g;
            for (itlj1g = rowFunctions.begin(); itlj1g != rowFunctions.end(); ++itlj1g) {
                LoopModel<Base>* loop = itlj1g->first;

                typename std::map<TapeVarType, std::map<size_t, GroupLoopRev2ColInfo<Base>* > >::const_iterator itj1g;
                for (itj1g = itlj1g->second.begin(); itj1g != itlj1g->second.end(); ++itj1g) {
                    const TapeVarType& jTape1 = itj1g->first;

                    typename std::map<size_t, GroupLoopRev2ColInfo<Base>* >::const_iterator itg;
                    for (itg = itj1g->second.begin(); itg != itj1g->second.end(); ++itg) {
                        size_t g = itg->first;
                        _cache << "         ";
                        generateFunctionNameLoopRev2(_cache, *loop, jTape1, g);
                        _cache << "(" << jrow << ", " << args << ");\n";
                    }
                }
            }

            /**
             * return all OK
             */
            _cache << "         return 0; // done\n";
        }
        _cache << "      default:\n"
                "         return 1; // error\n"
                "   };\n";

        _cache << "}\n";
        sources[functionRev2 + ".c"] = _cache.str();
        _cache.str("");

        /**
         * Sparsity
         */
        generateSparsity1DSource2(_name + "_" + FUNCTION_REVERSE_TWO_SPARSITY, elements);
        sources[_name + "_" + FUNCTION_REVERSE_TWO_SPARSITY + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateFunctionNameLoopRev2(std::ostringstream& cache,
                                                                     const LoopModel<Base>& loop,
                                                                     const TapeVarType& jTape1,
                                                                     size_t g) {
        cache << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO <<
                "_loop" << loop.getLoopId() <<
                "_indep" << jTape1.first << "_" << jTape1.second <<
                "_g" << g;
    }

}

#endif
