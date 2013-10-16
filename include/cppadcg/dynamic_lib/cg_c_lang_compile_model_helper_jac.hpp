#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_JAC_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_JAC_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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

    template<class Base>
    void CLangCompileModelHelper<Base>::generateJacobianSource(std::map<std::string, std::string>& sources) {
        const std::string jobName = "Jacobian";

        startingJob("operation graph for '" + jobName + "'");

        CodeHandler<Base> handler;
        handler.setJobTimer(this);

        vector<CGBase> indVars(_fun.Domain());
        handler.makeVariables(indVars);
        if (_x.size() > 0) {
            for (size_t i = 0; i < indVars.size(); i++) {
                indVars[i].setValue(_x[i]);
            }
        }

        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        vector<CGBase> jac(n * m);
        if (_jacMode == AUTOMATIC) {
            jac = _fun.Jacobian(indVars);
        } else if (_jacMode == FORWARD) {
            JacobianFor(_fun, indVars, jac);
        } else {
            JacobianRev(_fun, indVars, jac);
        }

        finishedJob();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_JACOBIAN);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("jac", "x", "var", "array"));

        handler.generateCode(code, langC, jac, *nameGen, _atomicFunctions, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseJacobianSource(std::map<std::string, std::string>& sources) throw (CGException) {
        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        /**
         * Determine the sparsity pattern
         */
        determineJacobianSparsity();

        bool forwardMode;

        if (_jacMode == AUTOMATIC) {
            if (_custom_jac.defined) {
                forwardMode = extra::estimateBestJacobianADMode(_jacSparsity.rows, _jacSparsity.cols);
            } else {
                forwardMode = n <= m;
            }
        } else {
            forwardMode = _jacMode == FORWARD;
        }

        /**
         * call the appropriate method for source code generation
         */
        if (_forwardOne && forwardMode) {
            generateSparseJacobianForRevSource(sources, true);
        } else if (_reverseOne && !forwardMode) {
            generateSparseJacobianForRevSource(sources, false);
        } else {
            generateSparseJacobianSource(sources, forwardMode);
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseJacobianSource(std::map<std::string, std::string>& sources,
                                                                     bool forward) throw (CGException) {
        const std::string jobName = "sparse Jacobian";

        //size_t m = _fun.Range();
        size_t n = _fun.Domain();

        startingJob("operation graph for '" + jobName + "'");

        CodeHandler<Base> handler;
        handler.setJobTimer(this);

        vector<CGBase> indVars(n);
        handler.makeVariables(indVars);
        if (_x.size() > 0) {
            for (size_t i = 0; i < n; i++) {
                indVars[i].setValue(_x[i]);
            }
        }

        vector<CGBase> jac(_jacSparsity.rows.size());
        if (_loopTapes.empty()) {
            //printSparsityPattern(_jacSparsity.sparsity, "jac sparsity");
            CppAD::sparse_jacobian_work work;
            if (forward) {
                _fun.SparseJacobianForward(indVars, _jacSparsity.sparsity, _jacSparsity.rows, _jacSparsity.cols, jac, work);
            } else {
                _fun.SparseJacobianReverse(indVars, _jacSparsity.sparsity, _jacSparsity.rows, _jacSparsity.cols, jac, work);
            }

        } else {
            jac = prepareSparseJacobianWithLoops(handler, indVars, forward);
        }

        finishedJob();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_SPARSE_JACOBIAN);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("jac", "x", "var", "array"));

        handler.generateCode(code, langC, jac, *nameGen, _atomicFunctions, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseJacobianForRevSource(std::map<std::string, std::string>& sources,
                                                                           bool forward) {
        //size_t m = _fun.Range();
        //size_t n = _fun.Domain();
        using namespace std;

        map<size_t, std::vector<size_t> > elements;
        map<size_t, std::vector<set<size_t> > > userJacElLocation; // maps each element to its position in the user jacobian
        string functionRevFor, revForSuffix;
        if (forward) {
            // elements[var]{equations}
            for (size_t e = 0; e < _jacSparsity.rows.size(); e++) {
                elements[_jacSparsity.cols[e]].push_back(_jacSparsity.rows[e]);
            }
            userJacElLocation = determineOrderByCol(elements, _jacSparsity);
            _cache.str("");
            _cache << _name << "_" << FUNCTION_SPARSE_FORWARD_ONE;
            functionRevFor = _cache.str();
            revForSuffix = "indep";
        } else {
            // elements[equation]{vars}
            for (size_t e = 0; e < _jacSparsity.rows.size(); e++) {
                elements[_jacSparsity.rows[e]].push_back(_jacSparsity.cols[e]);
            }
            userJacElLocation = determineOrderByRow(elements, _jacSparsity);
            _cache.str("");
            _cache << _name << "_" << FUNCTION_SPARSE_REVERSE_ONE;
            functionRevFor = _cache.str();
            revForSuffix = "dep";
        }


        /**
         * determine to which functions we can provide the jacobian row/column
         * directly without needing a temporary array (compressed)
         */
        map<size_t, bool> ordered;
        map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t index = it->first;
            const std::vector<size_t>& els = it->second;
            const std::vector<set<size_t> >& location = userJacElLocation.at(index);
            assert(els.size() == location.size());
            assert(els.size() > 0);

            bool passed = true;
            size_t jacArrayStart = *location[0].begin();
            for (size_t e = 0; e < els.size(); e++) {
                if (location[e].size() > 1) {
                    passed = false; // too many elements
                    break;
                }
                if (*location[e].begin() != jacArrayStart + e) {
                    passed = false; // wrong order
                    break;
                }
            }
            ordered[index] = passed;
        }
        assert(elements.size() == ordered.size());

        size_t maxCompressedSize = 0;
        map<size_t, bool>::const_iterator itOrd;
        for (it = elements.begin(), itOrd = ordered.begin(); it != elements.end(); ++it, ++itOrd) {
            if (it->second.size() > maxCompressedSize && !itOrd->second)
                maxCompressedSize = it->second.size();
        }

        if (!_loopTapes.empty()) {
            /**
             * with loops
             */
            if (forward) {
                generateSparseJacobianWithLoopsSourceFromForRev(sources, userJacElLocation, ordered, maxCompressedSize,
                                                                FUNCTION_SPARSE_FORWARD_ONE, "indep", "jcol",
                                                                _nonLoopFor1Elements, _loopFor1Groups,
                                                                generateFunctionNameLoopFor1);
            } else {
                generateSparseJacobianWithLoopsSourceFromForRev(sources, userJacElLocation, ordered, maxCompressedSize,
                                                                FUNCTION_SPARSE_REVERSE_ONE, "dep", "jrow",
                                                                _nonLoopRev1Elements, _loopRev1Groups,
                                                                generateFunctionNameLoopRev1);
            }
            return;
        }

        _cache.str("");
        _cache << _name << "_" << FUNCTION_SPARSE_JACOBIAN;
        string model_function(_cache.str());

        CLanguage<Base> langC(_baseTypeName);
        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                "\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n\n";
        generateFunctionDeclarationSource(_cache, functionRevFor, revForSuffix, elements, argsDcl);
        _cache << "\n"
                "void " << model_function << "(" << argsDcl << ") {\n"
                "   " << _baseTypeName << " const * inLocal[2];\n"
                "   " << _baseTypeName << " inLocal1 = 1;\n"
                "   " << _baseTypeName << " * outLocal[1];\n"
                "   " << _baseTypeName << " compressed[" << maxCompressedSize << "];\n"
                "   " << _baseTypeName << " * jac = out[0];\n"
                "\n"
                "   inLocal[0] = in[0];\n"
                "   inLocal[1] = &inLocal1;\n"
                "   outLocal[0] = compressed;";

        langC.setArgumentIn("inLocal");
        langC.setArgumentOut("outLocal");
        std::string argsLocal = langC.generateDefaultFunctionArguments();

        bool lastCompressed = true;
        for (it = elements.begin(), itOrd = ordered.begin(); it != elements.end(); ++it, ++itOrd) {
            size_t index = it->first;
            const std::vector<size_t>& els = it->second;
            const std::vector<set<size_t> >& location = userJacElLocation.at(index);
            assert(els.size() == location.size());

            _cache << "\n";
            if (itOrd->second) {
                _cache << "   outLocal[0] = &jac[" << *location[0].begin() << "];\n";
            } else if (!lastCompressed) {
                _cache << "   outLocal[0] = compressed;\n";
            }
            _cache << "   " << functionRevFor << "_" << revForSuffix << index << "(" << argsLocal << ");\n";
            if (!itOrd->second) {
                for (size_t e = 0; e < els.size(); e++) {
                    _cache << "   ";
                    set<size_t>::const_iterator itl;
                    for (itl = location[e].begin(); itl != location[e].end(); ++itl) {
                        _cache << "jac[" << (*itl) << "] = ";
                    }
                    _cache << "compressed[" << e << "];\n";
                }
            }
            lastCompressed = !itOrd->second;
        }

        _cache << "\n"
                "}\n";
        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::determineJacobianSparsity() {
        if (_jacSparsity.sparsity.size() > 0) {
            return;
        }

        /**
         * Determine the sparsity pattern
         */
        _jacSparsity.sparsity = extra::jacobianSparsitySet<SparsitySetType, CGBase> (_fun);

        if (!_custom_jac.defined) {
            extra::generateSparsityIndexes(_jacSparsity.sparsity, _jacSparsity.rows, _jacSparsity.cols);

        } else {
            _jacSparsity.rows = _custom_jac.row;
            _jacSparsity.cols = _custom_jac.col;
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateJacobianSparsitySource(std::map<std::string, std::string>& sources) {
        determineJacobianSparsity();

        generateSparsity2DSource(_name + "_" + FUNCTION_JACOBIAN_SPARSITY, _jacSparsity);
        sources[_name + "_" + FUNCTION_JACOBIAN_SPARSITY + ".c"] = _cache.str();
        _cache.str("");
    }

}

#endif
