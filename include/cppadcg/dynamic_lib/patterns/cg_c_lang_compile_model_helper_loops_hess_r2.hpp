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

    namespace loops {

        class RowGroup {
        public:
            std::auto_ptr<IndexPattern> jrowPattern;
            std::auto_ptr<IndexPattern> hessStartLocPattern;
        };
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseHessianWithLoopsSourceFromRev2(std::map<std::string, std::string>& sources,
                                                                                     const std::map<size_t, std::vector<std::set<size_t> > >& userHessElLocation,
                                                                                     const std::map<size_t, bool>& jrowOrdered,
                                                                                     size_t maxCompressedSize) {
        //size_t m = _fun->Range();
        //size_t n = _fun->Domain();
        using namespace std;
        using namespace CppAD::loops;

        /**
         * determine to which functions we can provide the hessian row directly
         * without needing a temporary array (compressed)
         */
        SmartVectorPointer<RowGroup> garbage;
        map<size_t, map<LoopModel<Base>*, map<size_t, RowGroup*> > > loopCalls;

        IndexDclrOperationNode<Base> indexIt("it");
        const string& itName = *indexIt.getName();

        std::vector<size_t> hessRowStart;
        /**
         * Determine jrow index patterns and
         * hessian row start patterns
         */
        std::vector<size_t> localit2jrows;
        typename map<LoopModel<Base>*, map<size_t, map<size_t, set<size_t> > > >::const_iterator itlge;
        for (itlge = _loopRev2Groups.begin(); itlge != _loopRev2Groups.end(); ++itlge) {
            LoopModel<Base>* loop = itlge->first;

            garbage.v.reserve(garbage.v.size() + itlge->second.size());

            map<size_t, map<size_t, set<size_t> > >::const_iterator itg;
            for (itg = itlge->second.begin(); itg != itlge->second.end(); ++itg) {
                size_t group = itg->first;
                const map<size_t, set<size_t> >& jrows2e = itg->second;

                // group by number of iterations
                std::auto_ptr<RowGroup> data(new RowGroup());

                /**
                 * jrow pattern
                 */
                mapKeys(jrows2e, localit2jrows);
                data->jrowPattern.reset(IndexPattern::detect(localit2jrows));

                /**
                 * hessian row start pattern
                 */
                bool ordered = true;
                for (size_t l = 0; l < localit2jrows.size(); l++) {
                    if (!jrowOrdered.at(localit2jrows[l])) {
                        ordered = false;
                        break;
                    }
                }

                if (ordered) {
                    hessRowStart.resize(localit2jrows.size());

                    for (size_t l = 0; l < localit2jrows.size(); l++) {
                        const std::vector<set<size_t> >& location = userHessElLocation.at(localit2jrows[l]);
                        hessRowStart[l] = *location[0].begin();
                    }

                    data->hessStartLocPattern.reset(IndexPattern::detect(hessRowStart));
                }

                // group by number of iterations
                loopCalls[localit2jrows.size()][loop][group] = data.get();
                garbage.v.push_back(data.release());
            }
        }

        string model_function = _name + "_" + FUNCTION_SPARSE_HESSIAN;
        string functionRev2 = _name + "_" + FUNCTION_SPARSE_REVERSE_TWO;
        string nlRev2Suffix = "noloop_indep";

        CLanguage<Base> langC(_baseTypeName);
        string loopFArgs = "inLocal, outLocal, " + langC.getArgumentAtomic();
        string argsDcl = langC.generateDefaultFunctionArgumentsDcl();


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
                "   unsigned long " << itName << ";\n"
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
        string argsLocal = langC.generateDefaultFunctionArguments();

        bool lastCompressed = true;
        map<size_t, set<size_t> >::const_iterator it;
        for (it = _nonLoopRev2Elements.begin(); it != _nonLoopRev2Elements.end(); ++it) {
            size_t index = it->first;
            const set<size_t>& elPos = it->second;
            const std::vector<set<size_t> >& location = userHessElLocation.at(index);
            assert(elPos.size() <= location.size()); // it can be lower because not all elements have to be assigned
            assert(elPos.size() > 0);
            bool rowOrdered = jrowOrdered.at(index);

            _cache << "\n";
            if (rowOrdered) {
                _cache << "   outLocal[0] = &hess[" << *location[0].begin() << "];\n";
            } else if (!lastCompressed) {
                _cache << "   outLocal[0] = compressed;\n";
            }
            _cache << "   " << functionRev2 << "_" << nlRev2Suffix << index << "(" << argsLocal << ");\n";
            if (!rowOrdered) {
                for (set<size_t>::const_iterator itEl = elPos.begin(); itEl != elPos.end(); ++itEl) {
                    size_t e = *itEl;
                    _cache << "   ";
                    set<size_t>::const_iterator itl;
                    for (itl = location[e].begin(); itl != location[e].end(); ++itl) {
                        _cache << "hess[" << (*itl) << "] += compressed[" << e << "];\n";
                    }
                }
            }
            lastCompressed = !rowOrdered;
        }

        /**
         * loop related values
         */
        typename map<size_t, map<LoopModel<Base>*, map<size_t, RowGroup*> > >::const_iterator itItlg;
        typename map<LoopModel<Base>*, map<size_t, RowGroup*> >::const_iterator itlg;
        typename map<size_t, RowGroup*>::const_iterator itg;

        lastCompressed = true;
        for (itItlg = loopCalls.begin(); itItlg != loopCalls.end(); ++itItlg) {
            size_t itCount = itItlg->first;
            if (itCount > 1) {
                _cache << "   for(" << itName << " = 0; " << itName << " < " << itCount << "; " << itName << "++) {\n";
            }

            for (itlg = itItlg->second.begin(); itlg != itItlg->second.end(); ++itlg) {
                LoopModel<Base>& loop = *itlg->first;

                for (itg = itlg->second.begin(); itg != itlg->second.end(); ++itg) {
                    size_t g = itg->first;
                    RowGroup* group = itg->second;

                    if (itCount > 1) {
                        _cache << "   "; //identation
                    }
                    if (group->hessStartLocPattern.get() != NULL) {
                        // determine hessRowStart = f(it)
                        _cache << "   outLocal[0] = &hess[" << CLanguage<Base>::indexPattern2String(*group->hessStartLocPattern, indexIt) << "];\n";
                    } else if (!lastCompressed) {
                        _cache << "   outLocal[0] = compressed;\n";
                    }

                    if (itCount > 1) {
                        _cache << "      jrow = " << CLanguage<Base>::indexPattern2String(*group->jrowPattern, indexIt) << ";\n";
                        _cache << "      ";
                        generateFunctionNameLoopRev2(_cache, loop, g);
                        _cache << "(jrow, " << loopFArgs << ");\n";
                    } else {
                        size_t jrow = _loopRev2Groups[&loop][g].begin()->first; // only one jrow
                        _cache << "   ";
                        generateFunctionNameLoopRev2(_cache, loop, g);
                        _cache << "(" << jrow << ", " << loopFArgs << ");\n";
                    }

                    if (group->hessStartLocPattern.get() == NULL) {
                        throw CGException("Not implemented yet");
                        _cache << "   outLocal[0] = compressed;\n"; //////////////////////
                    }

                    _cache << "\n";
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

        std::string argsDcl = langC.generateFunctionArgumentsDcl();
        std::string argsDclLoop = "unsigned long jrow, " + argsDcl;

        typename std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >::const_iterator itlg;
        for (itlg = _loopRev2Groups.begin(); itlg != _loopRev2Groups.end(); ++itlg) {
            const LoopModel<Base>& loop = *itlg->first;

            typename std::map<size_t, std::map<size_t, std::set<size_t> > >::const_iterator itg;
            for (itg = itlg->second.begin(); itg != itlg->second.end(); ++itg) {
                size_t group = itg->first;

                cache << "void ";
                generateFunctionNameLoopRev2(cache, loop, group);
                cache << "(" << argsDclLoop << ");\n";
            }
        }
    }

}

#endif