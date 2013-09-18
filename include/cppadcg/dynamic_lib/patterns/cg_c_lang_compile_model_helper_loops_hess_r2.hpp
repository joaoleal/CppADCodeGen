#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_HESS_R2_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_HESS_R2_INCLUDED
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
    void CLangCompileModelHelper<Base>::generateSparseHessianWithLoopsSourceFromRev2(std::map<std::string, std::string>& sources,
                                                                                     const std::map<size_t, std::vector<std::set<size_t> > >& userHessElLocation,
                                                                                     const std::map<size_t, bool>& jrowOrdered,
                                                                                     size_t maxCompressedSize) {
        //size_t m = _fun->Range();
        //size_t n = _fun->Domain();
        using namespace std;

        /**
         * determine to which functions we can provide the hessian row directly
         * without needing a temporary array (compressed)
         */

        std::map<size_t, std::map<LoopModel<Base>*, std::map<TapeVarType, map<size_t, GroupLoopRev2ColInfo<Base>*> > > > loopCalls;

        Index indexIt("it");
        const std::string& itName = indexIt.getName();

        std::vector<size_t> jrows, hessRowStart;
        /**
         * Determine jrow index patterns and
         * hessian row start patterns
         */
        typename std::map<LoopModel<Base>*, std::map<TapeVarType, vector<GroupLoopRev2ColInfo<Base>* > > >::const_iterator itljg;
        typename std::map<TapeVarType, vector<GroupLoopRev2ColInfo<Base>* > >::const_iterator itjg;

        for (itljg = _loopRev2Groups.begin(); itljg != _loopRev2Groups.end(); ++itljg) {
            LoopModel<Base>* loop = itljg->first;

            for (itjg = itljg->second.begin(); itjg != itljg->second.end(); ++itjg) {
                const TapeVarType& jTape1 = itjg->first;
                const vector<GroupLoopRev2ColInfo<Base>*>& groups = itjg->second;
                for (size_t g = 0; g < groups.size(); g++) {
                    GroupLoopRev2ColInfo<Base>* group = groups[g];
                    /**
                     * jrow pattern
                     */
                    mapKeys(group->jrow2localIt2ModelIt, jrows);

                    group->jrowPattern = IndexPattern::detect(indexIt, jrows);

                    /**
                     * hessian row start pattern
                     */
                    bool ordered = true;
                    std::map<size_t, std::map<size_t, size_t> >::const_iterator itJrow;
                    for (itJrow = group->jrow2localIt2ModelIt.begin(); itJrow != group->jrow2localIt2ModelIt.end(); ++itJrow) {
                        if (!jrowOrdered.at(itJrow->first)) {
                            ordered = false;
                            break;
                        }
                    }

                    if (ordered) {
                        hessRowStart.resize(jrows.size());

                        for (size_t e = 0; e < jrows.size(); e++) {
                            const std::vector<set<size_t> >& location = userHessElLocation.at(jrows[e]);
                            hessRowStart[e] = *location[0].begin();
                        }

                        group->hessStartLocPattern = IndexPattern::detect(indexIt, hessRowStart);
                    }

                    // group by number of iterations
                    loopCalls[jrows.size()][loop][jTape1][g] = group;
                }
            }
        }

        string model_function = _name + "_" + FUNCTION_SPARSE_HESSIAN;
        string functionRev2 = _name + "_" + FUNCTION_SPARSE_REVERSE_TWO;
        string nlRev2Suffix = "noloop_indep";

        CLanguage<Base> langC(_baseTypeName);
        std::string loopFArgs = "inLocal, outLocal, " + langC.getArgumentAtomic();
        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n\n";
        generateFunctionDeclarationSource(_cache, functionRev2, nlRev2Suffix, _nonLoopRev2Elements, argsDcl);
        generateFunctionDeclarationSourceLoopRev2(_cache, langC);

        _cache << "\n"
                "void " << model_function << "(" << argsDcl << ") {\n"
                "   " << _baseTypeName << " const * inLocal[3];\n"
                "   " << _baseTypeName << " inLocal1 = 1;\n"
                "   " << _baseTypeName << " * outLocal[1];\n"
                "   unsigned long " << indexIt.getName() << ";\n"
                "   unsigned long jrow;\n";
        if (maxCompressedSize > 0) {
            _cache << "   " << _baseTypeName << " compressed[" << maxCompressedSize << "];\n";
        }
        _cache << "   " << _baseTypeName << " * hess = out[0];\n"
                "\n"
                "   inLocal[0] = in[0];\n"
                "   inLocal[1] = &inLocal1;\n"
                "   inLocal[2] = in[1];\n";
        if (maxCompressedSize > 0) {
            _cache << "   outLocal[0] = compressed;\n";
        }
        _cache << "\n";

        /**
         * zero the hessian
         */
        _cache << "   for(" << itName << " = 0; " << itName << " < " << _hessSparsity.rows.size() << "; " << itName << "++) {\n"
                "      hess[" << itName << "] = 0;\n"
                "   }\n";

        /**
         * contributions from equations NOT belonging to loops
         * (must come before the loop related values because of the assigments)
         */
        langC.setArgumentIn("inLocal");
        langC.setArgumentOut("outLocal");
        std::string argsLocal = langC.generateDefaultFunctionArguments();

        bool lastCompressed = true;
        map<size_t, std::vector<Compressed2JColType> >::const_iterator it;
        for (it = _nonLoopRev2Elements.begin(); it != _nonLoopRev2Elements.end(); ++it) {
            size_t index = it->first;
            const std::vector<Compressed2JColType>& els = it->second;
            const std::vector<set<size_t> >& location = userHessElLocation.at(index);
            assert(els.size() <= location.size()); // it can be lower because not all elements have to be assigned
            assert(els.size() > 0);
            bool rowOrdered = jrowOrdered.at(index);

            _cache << "\n";
            if (rowOrdered) {
                _cache << "   outLocal[0] = &hess[" << *location[0].begin() << "];\n";
            } else if (!lastCompressed) {
                _cache << "   outLocal[0] = compressed;\n";
            }
            _cache << "   " << functionRev2 << "_" << nlRev2Suffix << index << "(" << argsLocal << ");\n";
            if (!rowOrdered) {
                for (size_t e = 0; e < els.size(); e++) {
                    _cache << "   ";
                    set<size_t>::const_iterator itl;
                    for (itl = location[e].begin(); itl != location[e].end(); ++itl) {
                        _cache << "hess[" << (*itl) << "] += compressed[" << e << "];\n";
                    }
                }
            }
            lastCompressed = !rowOrdered;
        }

        _cache << "\n";

        /**
         * loop related values
         */
        typename std::map<size_t, std::map<LoopModel<Base>*, std::map<TapeVarType, std::map<size_t, GroupLoopRev2ColInfo<Base>*> > > >::const_iterator itItljg2;
        typename std::map<LoopModel<Base>*, std::map<TapeVarType, std::map<size_t, GroupLoopRev2ColInfo<Base>*> > >::const_iterator itljg2;
        typename std::map<TapeVarType, std::map<size_t, GroupLoopRev2ColInfo<Base>*> >::const_iterator itjg2;
        typename std::map<size_t, GroupLoopRev2ColInfo<Base>*>::const_iterator itg;

        lastCompressed = true;
        for (itItljg2 = loopCalls.begin(); itItljg2 != loopCalls.end(); ++itItljg2) {
            size_t itCount = itItljg2->first;
            if (itCount > 1) {
                _cache << "   for(" << itName << " = 0; " << itName << " < " << itCount << "; " << itName << "++) {";
            }

            for (itljg2 = itItljg2->second.begin(); itljg2 != itItljg2->second.end(); ++itljg2) {
                const LoopModel<Base>& loop = *itljg2->first;
                for (itjg2 = itljg2->second.begin(); itjg2 != itljg2->second.end(); ++itjg2) {
                    const TapeVarType& jTape1 = itjg2->first;
                    for (itg = itjg2->second.begin(); itg != itjg2->second.end(); ++itg) {
                        size_t g = itg->first;
                        GroupLoopRev2ColInfo<Base>* group = itg->second;

                        _cache << "\n";
                        if (itCount > 1) {
                            _cache << "   "; //identation
                        }
                        if (group->hessStartLocPattern != NULL) {
                            // determine hessRowStart = f(it)
                            _cache << "   outLocal[0] = &hess[" << CLanguage<Base>::createIndexPattern(*group->hessStartLocPattern) << "];\n";
                        } else if (!lastCompressed) {
                            _cache << "   outLocal[0] = compressed;\n";
                        }

                        if (itCount > 1) {
                            _cache << "      jrow = " << CLanguage<Base>::createIndexPattern(*group->jrowPattern) << ";\n";
                            _cache << "      ";
                            generateFunctionNameLoopRev2(_cache, loop, jTape1, g);
                            _cache << "(jrow, " << loopFArgs << ");\n";
                        } else {
                            size_t jrow = group->jRow2iterations2Loc.begin()->first; // only one jrow
                            _cache << "   ";
                            generateFunctionNameLoopRev2(_cache, loop, jTape1, g);
                            _cache << "(" << jrow << ", " << loopFArgs << ");\n";
                        }

                        if (group->hessStartLocPattern == NULL) {
                            throw CGException("Not implemented yet");
                            _cache << "   outLocal[0] = compressed;\n"; //////////////////////
                        }
                    }
                }
            }

            if (itCount > 1) {
                _cache << "   }\n";
            }
        }

        _cache << "\n"
                "}\n";
        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateFunctionDeclarationSourceLoopRev2(std::ostringstream& cache,
                                                                                  CLanguage<Base>& langC) {

        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();
        std::string argsDclLoop = "unsigned long jrow, " + argsDcl;

        typename std::map<LoopModel<Base>*, std::map<TapeVarType, vector<GroupLoopRev2ColInfo<Base>* > > >::const_iterator itljg;
        for (itljg = _loopRev2Groups.begin(); itljg != _loopRev2Groups.end(); ++itljg) {
            const LoopModel<Base>& loop = *itljg->first;
            typename std::map<TapeVarType, vector<GroupLoopRev2ColInfo<Base>* > >::const_iterator itjg;
            for (itjg = itljg->second.begin(); itjg != itljg->second.end(); ++itjg) {
                const TapeVarType& jTape1 = itjg->first;
                const vector<GroupLoopRev2ColInfo<Base>*>& groups = itjg->second;
                for (size_t g = 0; g < groups.size(); g++) {
                    cache << "void ";
                    generateFunctionNameLoopRev2(cache, loop, jTape1, g);
                    cache << "(" << argsDclLoop << ");\n";
                }
            }
        }
    }

}

#endif