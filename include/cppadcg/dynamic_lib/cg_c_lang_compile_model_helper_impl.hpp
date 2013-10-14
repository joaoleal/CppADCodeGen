#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_IMPL_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_IMPL_INCLUDED
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

#include <typeinfo>
#include <memory>

namespace CppAD {

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_FORWAD_ZERO = "forward_zero";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_JACOBIAN = "jacobian";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_HESSIAN = "hessian";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_FORWARD_ONE = "forward_one";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_REVERSE_ONE = "reverse_one";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_REVERSE_TWO = "reverse_two";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_SPARSE_JACOBIAN = "sparse_jacobian";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_SPARSE_HESSIAN = "sparse_hessian";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_JACOBIAN_SPARSITY = "jacobian_sparsity";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_HESSIAN_SPARSITY = "hessian_sparsity";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_HESSIAN_SPARSITY2 = "hessian_sparsity2";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_SPARSE_FORWARD_ONE = "sparse_forward_one";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_SPARSE_REVERSE_ONE = "sparse_reverse_one";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_SPARSE_REVERSE_TWO = "sparse_reverse_two";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_FORWARD_ONE_SPARSITY = "forward_one_sparsity";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_REVERSE_ONE_SPARSITY = "reverse_one_sparsity";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_REVERSE_TWO_SPARSITY = "sparse_reverse_two_sparsity";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_INFO = "info";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_ATOMIC_FUNC_NAMES = "atomic_functions";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::CONST = "const";

    template<class Base>
    VariableNameGenerator<Base>* CLangCompileModelHelper<Base>::createVariableNameGenerator(const std::string& depName,
                                                                                            const std::string& indepName,
                                                                                            const std::string& tmpName,
                                                                                            const std::string& tmpArrayName) {
        return new CLangDefaultVariableNameGenerator<Base > (depName, indepName, tmpName, tmpArrayName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::compileSources(CLangCompiler<Base>& compiler,
                                                       bool posIndepCode) throw (CGException) {
        generateLoops();

        std::map<std::string, std::string> sources;
        if (_zero) {
            generateZeroSource(sources);
        }

        if (_jacobian) {
            generateJacobianSource(sources);
        }

        if (_hessian) {
            generateHessianSource(sources);
        }

        if (_forwardOne) {
            generateSparseForwardOneSources(sources);
            generateForwardOneSources(sources);
        }

        if (_reverseOne) {
            generateSparseReverseOneSources(sources);
            generateReverseOneSources(sources);
        }

        if (_reverseTwo) {
            generateSparseReverseTwoSources(sources);
            generateReverseTwoSources(sources);
        }

        if (_sparseJacobian) {
            generateSparseJacobianSource(sources);
        }

        if (_sparseHessian) {
            generateSparseHessianSource(sources);
        }

        if (_sparseJacobian || _forwardOne || _reverseOne) {
            generateJacobianSparsitySource(sources);
        }

        if (_sparseHessian || _reverseTwo) {
            generateHessianSparsitySource(sources);
        }

        generateInfoSource(sources);

        generateAtomicFuncNames(sources);

        compiler.compileSources(sources, posIndepCode, true);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateLoops() throw (CGException) {
        if (_relatedDepCandidates.empty()) {
            return; //nothing to do
        }

        startingGraphCreation("Loop detection");

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        std::vector<CGBase> xx(_fun.Domain());
        handler.makeVariables(xx);
        if (_x.size() > 0) {
            for (size_t i = 0; i < xx.size(); i++) {
                xx[i].setValue(_x[i]);
            }
        }

        std::vector<CGBase> yy = _fun.Forward(0, xx);

        DependentPatternMatcher<Base> matcher(_relatedDepCandidates, yy, xx);
        matcher.generateTapes(_funNoLoops, _loopTapes);

        if (_verbose) {
            std::cout << "equation patterns: " << matcher.getEquationPatterns().size() << std::endl;
            std::cout << "loops: " << matcher.getLoops().size() << std::endl;
        }

        finishedGraphCreation();
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateInfoSource(std::map<std::string, std::string>& sources) {
        const char* localBaseName = typeid (Base).name();

        std::string funcName = _name + "_" + FUNCTION_INFO;

        std::auto_ptr<VariableNameGenerator< Base > > nameGen(createVariableNameGenerator("y", "x", "var", "array"));

        _cache.str("");
        _cache << "void " << funcName << "(const char** baseName, unsigned long* m, unsigned long* n, unsigned int* indCount, unsigned int* depCount) {\n"
                "   *baseName = \"" << _baseTypeName << "  " << localBaseName << "\";\n"
                "   *m = " << _fun.Range() << ";\n"
                "   *n = " << _fun.Domain() << ";\n"
                "   *depCount = " << nameGen->getDependent().size() << "; // number of dependent array variables\n"
                "   *indCount = " << nameGen->getIndependent().size() << "; // number of independent array variables\n"
                "}\n\n";

        sources[funcName + ".c"] = _cache.str();
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateAtomicFuncNames(std::map<std::string, std::string>& sources) {
        std::string funcName = _name + "_" + FUNCTION_ATOMIC_FUNC_NAMES;
        size_t n = _atomicFunctions.size();
        _cache.str("");
        _cache << "void " << funcName << "(const char*** names, unsigned long* n) {\n"
                "   static const char* atomic[" << n << "] = {";
        for (size_t i = 0; i < n; i++) {
            if (i > 0)_cache << ", ";
            _cache << "\"" << _atomicFunctions[i] << "\"";
        }
        _cache << "};\n"
                "   *names = atomic;\n"
                "   *n = " << n << ";\n"
                "}\n\n";

        sources[funcName + ".c"] = _cache.str();
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateZeroSource(std::map<std::string, std::string>& sources) {
        const std::string jobName = "model (zero-order forward)";

        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        vector<CGBase> indVars(_fun.Domain());
        handler.makeVariables(indVars);
        if (_x.size() > 0) {
            for (size_t i = 0; i < indVars.size(); i++) {
                indVars[i].setValue(_x[i]);
            }
        }

        vector<CGBase> dep;

        if (_loopTapes.empty()) {
            dep = _fun.Forward(0, indVars);
        } else {
            /**
             * Contains loops
             */
            dep = prepareForward0WithLoops(handler, indVars);
        }

        finishedGraphCreation();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_FORWAD_ZERO);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("y", "x", "var", "array"));

        handler.generateCode(code, langC, dep, *nameGen, _atomicFunctions, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateJacobianSource(std::map<std::string, std::string>& sources) {
        const std::string jobName = "Jacobian";

        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

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

        finishedGraphCreation();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_JACOBIAN);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("jac", "x", "var", "array"));

        handler.generateCode(code, langC, jac, *nameGen, _atomicFunctions, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateHessianSource(std::map<std::string, std::string>& sources) {
        const std::string jobName = "Hessian";

        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        size_t m = _fun.Range();
        size_t n = _fun.Domain();


        // independent variables
        vector<CGBase> indVars(n);
        handler.makeVariables(indVars);
        if (_x.size() > 0) {
            for (size_t i = 0; i < n; i++) {
                indVars[i].setValue(_x[i]);
            }
        }

        // multipliers
        vector<CGBase> w(m);
        handler.makeVariables(w);
        if (_x.size() > 0) {
            for (size_t i = 0; i < m; i++) {
                w[i].setValue(Base(1.0));
            }
        }

        vector<CGBase> hess = _fun.Hessian(indVars, w);

        // make use of the symmetry of the Hessian in order to reduce operations
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < i; j++) {
                hess[i * n + j] = hess[j * n + i];
            }
        }

        finishedGraphCreation();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_HESSIAN);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("hess", "x", "var", "array"));
        CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), n);

        handler.generateCode(code, langC, hess, nameGenHess, _atomicFunctions, jobName);
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

        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

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

        finishedGraphCreation();

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

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseHessianSource(std::map<std::string, std::string>& sources) throw (CGException) {
        /**
         * Determine the sparsity pattern p for Hessian of w^T F
         */
        determineHessianSparsity();

        if (_reverseTwo) {
            generateSparseHessianSourceFromRev2(sources);
        } else {
            generateSparseHessianSourceDirectly(sources);
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseHessianSourceDirectly(std::map<std::string, std::string>& sources) throw (CGException) {
        const std::string jobName = "sparse Hessian";
        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        /**
         * we might have to consider a slightly different order than the one
         * specified by the user according to the available elements in the sparsity
         */
        std::vector<size_t> evalRows, evalCols;
        determineSecondOrderElements4Eval(evalRows, evalCols);

        std::map<size_t, std::map<size_t, size_t> > locations;
        for (size_t e = 0; e < evalRows.size(); e++) {
            size_t j1 = evalRows[e];
            size_t j2 = evalCols[e];
            std::map<size_t, std::map<size_t, size_t> >::iterator itJ1 = locations.find(j1);
            if (itJ1 == locations.end()) {
                locations[j1][j2] = e;
            } else {
                std::map<size_t, size_t>& j22e = itJ1->second;
                if (j22e.find(j2) == j22e.end()) {
                    j22e[j2] = e; // OK
                } else {
                    // repeated elements not allowed
                    std::ostringstream ss;
                    ss << "Repeated Hessian element requested: " << j1 << " " << j2;
                    throw CGException(ss.str());
                }
            }
        }

        // make use of the symmetry of the Hessian in order to reduce operations
        std::vector<size_t> lowerHessRows, lowerHessCols, lowerHessOrder;
        lowerHessRows.reserve(_hessSparsity.rows.size() / 2);
        lowerHessCols.reserve(lowerHessRows.size());
        lowerHessOrder.reserve(lowerHessRows.size());

        std::map<size_t, size_t> duplicates; // the elements determined using symmetry
        std::map<size_t, std::map<size_t, size_t> >::const_iterator itJ;
        std::map<size_t, size_t>::const_iterator itI;
        for (size_t e = 0; e < evalRows.size(); e++) {
            bool add = true;
            size_t i = evalRows[e];
            size_t j = evalCols[e];
            if (i < j) {
                // find the symmetric value
                itJ = locations.find(j);
                if (itJ != locations.end()) {
                    itI = itJ->second.find(i);
                    if (itI != itJ->second.end()) {
                        size_t eSim = itI->second;
                        duplicates[e] = eSim;
                        add = false; // symmetric value being determined
                    }
                }
            }

            if (add) {
                lowerHessRows.push_back(i);
                lowerHessCols.push_back(j);
                lowerHessOrder.push_back(e);
            }
        }

        /**
         * 
         */
        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        // independent variables
        vector<CGBase> indVars(n);
        handler.makeVariables(indVars);
        if (_x.size() > 0) {
            for (size_t i = 0; i < n; i++) {
                indVars[i].setValue(_x[i]);
            }
        }

        // multipliers
        vector<CGBase> w(m);
        handler.makeVariables(w);
        if (_x.size() > 0) {
            for (size_t i = 0; i < m; i++) {
                w[i].setValue(Base(1.0));
            }
        }

        vector<CGBase> hess(_hessSparsity.rows.size());
        if (_loopTapes.empty()) {
            CppAD::sparse_hessian_work work;
            vector<CGBase> lowerHess(lowerHessRows.size());
            _fun.SparseHessian(indVars, w, _hessSparsity.sparsity, lowerHessRows, lowerHessCols, lowerHess, work);

            for (size_t i = 0; i < lowerHessOrder.size(); i++) {
                hess[lowerHessOrder[i]] = lowerHess[i];
            }

            // make use of the symmetry of the Hessian in order to reduce operations
            std::map<size_t, size_t>::const_iterator it2;
            for (it2 = duplicates.begin(); it2 != duplicates.end(); ++it2) {
                hess[it2->first] = hess[it2->second];
            }
        } else {
            /**
             * with loops
             */
            hess = prepareSparseHessianWithLoops(handler, indVars, w,
                                                 lowerHessRows, lowerHessCols, lowerHessOrder,
                                                 duplicates);
        }

        finishedGraphCreation();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_SPARSE_HESSIAN);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("hess", "x", "var", "array"));
        CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), n);

        handler.generateCode(code, langC, hess, nameGenHess, _atomicFunctions, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseHessianSourceFromRev2(std::map<std::string, std::string>& sources) throw (CGException) {
        using namespace std;

        /**
         * we might have to consider a slightly different order than the one
         * specified by the user according to the available elements in the sparsity
         */
        std::vector<size_t> evalRows, evalCols;
        determineSecondOrderElements4Eval(evalRows, evalCols);

        // elements[var]{var}
        map<size_t, std::vector<size_t> > elements;
        for (size_t e = 0; e < evalRows.size(); e++) {
            elements[evalRows[e]].push_back(evalCols[e]);
        }

        // maps each element to its position in the user hessian
        map<size_t, std::vector<set<size_t> > > userHessElLocation = determineOrderByRow(elements, evalRows, evalCols);

        /**
         * determine to which functions we can provide the hessian row directly
         * without needing a temporary array (compressed)
         */
        map<size_t, bool> ordered;
        map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t index = it->first;
            const std::vector<size_t>& els = it->second;
            const std::vector<set<size_t> >& location = userHessElLocation.at(index);
            assert(els.size() == location.size());
            assert(els.size() > 0);

            bool passed = true;
            size_t hessRowStart = *location[0].begin();
            for (size_t e = 0; e < els.size(); e++) {
                if (location[e].size() > 1) {
                    passed = false; // too many elements
                    break;
                }
                if (*location[e].begin() != hessRowStart + e) {
                    passed = false; // wrong order
                    break;
                }
            }
            ordered[index] = passed;
        }
        assert(elements.size() == ordered.size());

        /**
         * determine the maximum size of the temporary array
         */
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
            generateSparseHessianWithLoopsSourceFromRev2(sources, userHessElLocation, ordered, maxCompressedSize);
            return;
        }

        string model_function = _name + "_" + FUNCTION_SPARSE_HESSIAN;
        string functionRev2 = _name + "_" + FUNCTION_SPARSE_REVERSE_TWO;
        string rev2Suffix = "indep";

        CLanguage<Base> langC(_baseTypeName);
        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n\n";
        generateFunctionDeclarationSource(_cache, functionRev2, rev2Suffix, elements, argsDcl);
        _cache << "\n"
                "void " << model_function << "(" << argsDcl << ") {\n"
                "   " << _baseTypeName << " const * inLocal[3];\n"
                "   " << _baseTypeName << " inLocal1 = 1;\n"
                "   " << _baseTypeName << " * outLocal[1];\n";
        if (maxCompressedSize > 0) {
            _cache << "   " << _baseTypeName << " compressed[" << maxCompressedSize << "];\n";
        }
        _cache << "   " << _baseTypeName << " * hess = out[0];\n"
                "\n"
                "   inLocal[0] = in[0];\n"
                "   inLocal[1] = &inLocal1;\n"
                "   inLocal[2] = in[1];\n";
        if (maxCompressedSize > 0) {
            _cache << "   outLocal[0] = compressed;";
        }

        langC.setArgumentIn("inLocal");
        langC.setArgumentOut("outLocal");
        std::string argsLocal = langC.generateDefaultFunctionArguments();
        bool lastCompressed = true;
        for (it = elements.begin(), itOrd = ordered.begin(); it != elements.end(); ++it, ++itOrd) {
            size_t index = it->first;
            const std::vector<size_t>& els = it->second;
            const std::vector<set<size_t> >& location = userHessElLocation.at(index);
            assert(els.size() == location.size());
            assert(els.size() > 0);

            _cache << "\n";
            if (itOrd->second) {
                _cache << "   outLocal[0] = &hess[" << *location[0].begin() << "];\n";
            } else if (!lastCompressed) {
                _cache << "   outLocal[0] = compressed;\n";
            }
            _cache << "   " << functionRev2 << "_" << rev2Suffix << index << "(" << argsLocal << ");\n";
            if (!itOrd->second) {
                for (size_t e = 0; e < els.size(); e++) {
                    _cache << "   ";
                    set<size_t>::const_iterator itl;
                    for (itl = location[e].begin(); itl != location[e].end(); ++itl) {
                        _cache << "hess[" << (*itl) << "] = ";
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
    void CLangCompileModelHelper<Base>::determineSecondOrderElements4Eval(std::vector<size_t>& evalRows,
                                                                          std::vector<size_t>& evalCols) {
        /**
         * Atomic functions migth not have all the elements and thus there may 
         * be no symmetry. This will explore symmetry in order to provide the
         * second order elements requested by the user.
         */
        evalRows.reserve(_hessSparsity.rows.size());
        evalCols.reserve(_hessSparsity.cols.size());

        for (size_t e = 0; e < _hessSparsity.rows.size(); e++) {
            size_t i = _hessSparsity.rows[e];
            size_t j = _hessSparsity.cols[e];
            if (_hessSparsity.sparsity[i].find(j) == _hessSparsity.sparsity[i].end() &&
                    _hessSparsity.sparsity[j].find(i) != _hessSparsity.sparsity[j].end()) {
                // only the symmetric value is available
                // (it can be caused by atomic functions which may only be providing a partial hessian)
                evalRows.push_back(j);
                evalCols.push_back(i);
            } else {
                evalRows.push_back(i);
                evalCols.push_back(j);
            }
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::determineHessianSparsity() {
        if (_hessSparsity.sparsity.size() > 0) {
            return;
        }

        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        /**
         * sparsity for the sum of the hessians of all equations
         */
        SparsitySetType r(n); // identity matrix
        for (size_t j = 0; j < n; j++)
            r[j].insert(j);
        SparsitySetType jac = _fun.ForSparseJac(n, r);

        SparsitySetType s(1);
        for (size_t i = 0; i < m; i++) {
            s[0].insert(i);
        }
        _hessSparsity.sparsity = _fun.RevSparseHes(n, s, false);
        //printSparsityPattern(_hessSparsity.sparsity, "hessian");

        if (_hessianByEquation || _reverseTwo) {
            /**
             * sparsity for the hessian of each equations
             */

            std::set<size_t>::const_iterator it;

            std::set<size_t> customVarsInHess;
            if (_custom_hess.defined) {
                customVarsInHess.insert(_custom_hess.row.begin(), _custom_hess.row.end());
                customVarsInHess.insert(_custom_hess.col.begin(), _custom_hess.col.end());

                r = SparsitySetType(n); //clear r
                for (it = customVarsInHess.begin(); it != customVarsInHess.end(); ++it) {
                    size_t j = *it;
                    r[j].insert(j);
                }
                jac = _fun.ForSparseJac(n, r);
            }

            /**
             * Coloring
             */
            const CppAD::vector<Color> colors = colorByRow(customVarsInHess, jac);

            /**
             * For each individual equation
             */
            _hessSparsities.resize(m);
            for (size_t i = 0; i < m; i++) {
                _hessSparsities[i].sparsity.resize(n);
            }

            for (size_t c = 0; c < colors.size(); c++) {
                const Color& color = colors[c];

                // first-order
                r = SparsitySetType(n); //clear r
                for (it = color.forbiddenRows.begin(); it != color.forbiddenRows.end(); ++it) {
                    size_t j = *it;
                    r[j].insert(j);
                }
                _fun.ForSparseJac(n, r);

                // second-order
                s[0].clear();
                const std::set<size_t>& equations = color.rows;
                std::set<size_t>::const_iterator itEq;
                for (itEq = equations.begin(); itEq != equations.end(); ++itEq) {
                    size_t i = *itEq;
                    s[0].insert(i);
                }

                SparsitySetType sparsityc = _fun.RevSparseHes(n, s, false);

                /**
                 * Retrieve the individual hessians for each equation
                 */
                const std::map<size_t, size_t>& var2Eq = color.column2Row;
                const std::set<size_t>& forbidden_c = color.forbiddenRows; //used variables
                std::set<size_t>::const_iterator itvar;
                for (itvar = forbidden_c.begin(); itvar != forbidden_c.end(); ++itvar) {
                    size_t j = *itvar;
                    if (sparsityc[j].size() > 0) {
                        size_t i = var2Eq.at(j);
                        _hessSparsities[i].sparsity[j].insert(sparsityc[j].begin(),
                                                              sparsityc[j].end());
                    }
                }

            }

            for (size_t i = 0; i < m; i++) {
                LocalSparsityInfo& hessSparsitiesi = _hessSparsities[i];

                if (!_custom_hess.defined) {
                    extra::generateSparsityIndexes(hessSparsitiesi.sparsity,
                                                   hessSparsitiesi.rows, hessSparsitiesi.cols);

                } else {
                    size_t nnz = _custom_hess.row.size();
                    for (size_t e = 0; e < nnz; e++) {
                        size_t i1 = _custom_hess.row[e];
                        size_t i2 = _custom_hess.col[e];
                        if (hessSparsitiesi.sparsity[i1].find(i2) != hessSparsitiesi.sparsity[i1].end()) {
                            hessSparsitiesi.rows.push_back(i1);
                            hessSparsitiesi.cols.push_back(i2);
                        }
                    }
                }
            }

        }

        if (!_custom_hess.defined) {
            extra::generateSparsityIndexes(_hessSparsity.sparsity,
                                           _hessSparsity.rows, _hessSparsity.cols);

        } else {
            _hessSparsity.rows = _custom_hess.row;
            _hessSparsity.cols = _custom_hess.col;
        }
    }

    template<class Base>
    vector<typename CLangCompileModelHelper<Base>::Color> CLangCompileModelHelper<Base>::colorByRow(const std::set<size_t>& columns,
                                                                                                    const SparsitySetType& sparsity) {
        CppAD::vector<Color> colors(sparsity.size()); // reserve the maximum size to avoid reallocating more space later

        std::set<size_t>::const_iterator it;

        /**
         * try not match the columns of each row to a color which did not have
         * those columns yet 
         */
        size_t c_used = 0;
        for (size_t i = 0; i < sparsity.size(); i++) {
            const std::set<size_t>& row = sparsity[i];
            if (row.size() == 0) {
                continue; //nothing to do
            }

            // consider only the columns present in the sparsity pattern
            std::set<size_t> rowReduced;
            if (_custom_hess.defined) {
                for (it = row.begin(); it != row.end(); ++it) {
                    size_t j = *it;
                    if (columns.find(j) != columns.end())
                        rowReduced.insert(j);
                }
            } else {
                rowReduced = row;
            }

            bool newColor = true;
            size_t colori;
            for (size_t c = 0; c < c_used; c++) {
                std::set<size_t>& forbidden_c = colors[c].forbiddenRows;
                if (!intersects(forbidden_c, rowReduced)) {
                    // no intersection
                    colori = c;
                    newColor = false;
                    forbidden_c.insert(rowReduced.begin(), rowReduced.end());
                    break;
                }
            }

            if (newColor) {
                colori = c_used;
                colors[c_used].forbiddenRows = rowReduced;
                c_used++;
            }

            colors[colori].rows.insert(i);

            for (it = rowReduced.begin(); it != rowReduced.end(); ++it) {
                size_t j = *it;
                colors[colori].column2Row[j] = i;
                colors[colori].row2Columns[i].insert(j);
            }
        }

        colors.resize(c_used); //reduce size
        return colors;
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateHessianSparsitySource(std::map<std::string, std::string>& sources) {
        determineHessianSparsity();

        generateSparsity2DSource(_name + "_" + FUNCTION_HESSIAN_SPARSITY, _hessSparsity);
        sources[_name + "_" + FUNCTION_HESSIAN_SPARSITY + ".c"] = _cache.str();
        _cache.str("");

        if (_hessianByEquation || _reverseTwo) {
            generateSparsity2DSource2(_name + "_" + FUNCTION_HESSIAN_SPARSITY2, _hessSparsities);
            sources[_name + "_" + FUNCTION_HESSIAN_SPARSITY2 + ".c"] = _cache.str();
            _cache.str("");
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseForwardOneSources(std::map<std::string, std::string>& sources) {
#ifndef NDEBUG
        size_t m = _fun.Range();
#endif
        size_t n = _fun.Domain();

        determineJacobianSparsity();

        // elements[var]{equations}
        std::map<size_t, std::vector<size_t> > elements;
        for (size_t e = 0; e < _jacSparsity.rows.size(); e++) {
            elements[_jacSparsity.cols[e]].push_back(_jacSparsity.rows[e]);
        }

        if (!_loopTapes.empty()) {
            /**
             * with loops
             */
            prepareSparseForwardOneWithLoops(sources, elements);
            return;
        }

        /**
         * Generate one function for each dependent variable
         */
        vector<CGBase> dxv(n);

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t j = it->first;
            const std::vector<size_t>& rows = it->second;

            _cache.str("");
            _cache << "model (forward one)";
            const std::string jobName = _cache.str();

            startingGraphCreation(jobName);

            CodeHandler<Base> handler;
            handler.setVerbose(_verbose);

            vector<CGBase> indVars(n);
            handler.makeVariables(indVars);
            if (_x.size() > 0) {
                for (size_t i = 0; i < n; i++) {
                    indVars[i].setValue(_x[i]);
                }
            }

            CGBase dx;
            handler.makeVariable(dx);
            if (_x.size() > 0) {
                dx.setValue(Base(1.0));
            }

            // TODO: consider caching the zero order coefficients somehow between calls
            _fun.Forward(0, indVars);
            dxv[j] = dx;
            vector<CGBase> dy = _fun.Forward(1, dxv);
            dxv[j] = Base(0);
            assert(dy.size() == m);

            vector<CGBase> dyCustom;
            std::vector<size_t>::const_iterator it2;
            for (it2 = rows.begin(); it2 != rows.end(); ++it2) {
                dyCustom.push_back(dy[*it2]);
            }

            finishedGraphCreation();

            CLanguage<Base> langC(_baseTypeName);
            langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
            _cache.str("");
            _cache << _name << "_" << FUNCTION_SPARSE_FORWARD_ONE << "_indep" << j;
            langC.setGenerateFunction(_cache.str());

            std::ostringstream code;
            std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("dy", "x", "var", "array"));
            CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), "dx", n);

            handler.generateCode(code, langC, dyCustom, nameGenHess, _atomicFunctions, jobName);
        }

        _cache.str("");

        generateGlobalDirectionalFunctionSource(FUNCTION_SPARSE_FORWARD_ONE,
                                                "indep",
                                                FUNCTION_FORWARD_ONE_SPARSITY,
                                                elements,
                                                sources);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateForwardOneSources(std::map<std::string, std::string>& sources) {

        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        _cache.str("");
        _cache << _name << "_" << FUNCTION_FORWARD_ONE;
        std::string model_function(_cache.str());

        CLanguage<Base> langC(_baseTypeName);
        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();
        std::string args = langC.generateDefaultFunctionArguments();

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n"
                "\n"
                "void " << _name << "_" << FUNCTION_SPARSE_FORWARD_ONE << "(unsigned long pos, " << argsDcl << ");\n"
                "void " << _name << "_" << FUNCTION_FORWARD_ONE_SPARSITY << "(unsigned long pos, unsigned long const** elements, unsigned long* nnz);\n"
                "\n"
                "int " << model_function << "("
                << _baseTypeName << " const tx[], "
                << _baseTypeName << " ty[], "
                << langC.generateArgumentAtomicDcl() << ") {\n"
                "   unsigned long ePos, ej, i, j, nnz, nnzMax;\n"
                "   unsigned long const* pos;\n"
                "   unsigned long* tyPos;\n"
                "   unsigned long nnzTy;\n"
                "   " << _baseTypeName << " const * in[2];\n"
                "   " << _baseTypeName << "* out[1];\n"
                "   " << _baseTypeName << " x[" << n << "];\n"
                "   " << _baseTypeName << "* compressed;\n"
                "\n"
                "   tyPos = 0;\n"
                "   nnzTy = 0;\n"
                "   nnzMax = 0;\n"
                "   for (j = 0; j < " << n << "; j++) {\n"
                "      if (tx[j * 2 + 1] != 0.0) {\n"
                "         nnzTy++;\n"
                "         tyPos = (unsigned long*) realloc(tyPos, nnzTy * sizeof(unsigned long));\n"
                "         tyPos[nnzTy - 1] = j;\n"
                "         " << _name << "_" << FUNCTION_FORWARD_ONE_SPARSITY << "(j, &pos, &nnz);\n"
                "         if(nnz > nnzMax)\n"
                "            nnzMax = nnz;\n"
                "      }\n"
                "   }\n"
                "   for (i = 0; i < " << m << "; i++) {\n"
                "      ty[i * 2 + 1] = 0;\n"
                "   }\n"
                "\n"
                "   if (nnzTy == 0) {\n"
                "      free(tyPos);\n"
                "      return 0; //nothing to do\n"
                "   }\n"
                "\n"
                "   compressed = (" << _baseTypeName << "*) malloc(nnzMax * sizeof(" << _baseTypeName << "));\n"
                "\n"
                "   for (j = 0; j < " << n << "; j++)\n"
                "      x[j] = tx[j * 2];\n"
                "\n"
                "   for (ej = 0; ej < nnzTy; ej++) {\n"
                "      j = tyPos[ej];\n"
                "      " << _name << "_" << FUNCTION_FORWARD_ONE_SPARSITY << "(j, &pos, &nnz);\n"
                "\n"
                "      in[0] = x;\n"
                "      in[1] = &tx[j * 2 + 1];\n"
                "      out[0] = compressed;\n"
                "      " << _name << "_" << FUNCTION_SPARSE_FORWARD_ONE << "(j, " << args << ");\n"
                "\n"
                "      for (ePos = 0; ePos < nnz; ePos++) {\n"
                "         ty[pos[ePos] * 2 + 1] += compressed[ePos];\n"
                "      }\n"
                "\n"
                "   }\n"
                "   free(compressed);\n"
                "   free(tyPos);\n"
                "   return 0;\n"
                "}\n";
        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseReverseOneSources(std::map<std::string, std::string>& sources) {
        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        determineJacobianSparsity();

        // elements[equation]{vars}
        std::map<size_t, std::vector<size_t> > elements;
        for (size_t e = 0; e < _jacSparsity.rows.size(); e++) {
            elements[_jacSparsity.rows[e]].push_back(_jacSparsity.cols[e]);
        }

        if (!_loopTapes.empty()) {
            /**
             * with loops
             */
            prepareSparseReverseOneWithLoops(sources, elements);
            return;
        }

        vector<CGBase> w(m);

        /**
         * Generate one function for each dependent variable
         */
        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t i = it->first;
            const std::vector<size_t>& cols = it->second;

            _cache.str("");
            _cache << "model (reverse one, dep " << i << ")";
            const std::string jobName = _cache.str();

            startingGraphCreation(jobName);

            CodeHandler<Base> handler;
            handler.setVerbose(_verbose);

            vector<CGBase> indVars(_fun.Domain());
            handler.makeVariables(indVars);
            if (_x.size() > 0) {
                for (size_t i = 0; i < n; i++) {
                    indVars[i].setValue(_x[i]);
                }
            }

            CGBase py;
            handler.makeVariable(py);
            if (_x.size() > 0) {
                py.setValue(Base(1.0));
            }

            // TODO: consider caching the zero order coefficients somehow between calls
            _fun.Forward(0, indVars);

            w[i] = py;
            vector<CGBase> dw = _fun.Reverse(1, w);
            assert(dw.size() == n);
            w[i] = Base(0);

            vector<CGBase> dwCustom;
            std::vector<size_t>::const_iterator it2;
            for (it2 = cols.begin(); it2 != cols.end(); ++it2) {
                dwCustom.push_back(dw[*it2]);
            }

            finishedGraphCreation();

            CLanguage<Base> langC(_baseTypeName);
            langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
            _cache.str("");
            _cache << _name << "_" << FUNCTION_SPARSE_REVERSE_ONE << "_dep" << i;
            langC.setGenerateFunction(_cache.str());

            std::ostringstream code;
            std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("dw", "x", "var", "array"));
            CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), "py", n);

            handler.generateCode(code, langC, dwCustom, nameGenHess, _atomicFunctions, jobName);
        }

        _cache.str("");

        generateGlobalDirectionalFunctionSource(FUNCTION_SPARSE_REVERSE_ONE,
                                                "dep",
                                                FUNCTION_REVERSE_ONE_SPARSITY,
                                                elements,
                                                sources);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateReverseOneSources(std::map<std::string, std::string>& sources) {
        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        _cache.str("");
        _cache << _name << "_" << FUNCTION_REVERSE_ONE;
        std::string model_function(_cache.str());
        _cache.str("");

        CLanguage<Base> langC(_baseTypeName);
        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();
        std::string args = langC.generateDefaultFunctionArguments();

        _cache << "#include <stdlib.h>\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n"
                "\n"
                "void " << _name << "_" << FUNCTION_SPARSE_REVERSE_ONE << "(unsigned long pos, " << argsDcl << ");\n"
                "void " << _name << "_" << FUNCTION_REVERSE_ONE_SPARSITY << "(unsigned long pos, unsigned long const** elements, unsigned long* nnz);\n"
                "\n"
                "int " << model_function << "("
                << _baseTypeName << " const x[], "
                << _baseTypeName << " const ty[],"
                << _baseTypeName << "  px[], "
                << _baseTypeName << " const py[], "
                << langC.generateArgumentAtomicDcl() << ") {\n"
                "   unsigned long ei, ePos, i, j, nnz, nnzMax;\n"
                "   unsigned long const* pos;\n"
                "   unsigned long* pyPos;\n"
                "   unsigned long nnzPy;\n"
                "   " << _baseTypeName << " const * in[2];\n"
                "   " << _baseTypeName << "* out[1];\n"
                "   " << _baseTypeName << "* compressed;\n"
                "\n"
                "   pyPos = 0;\n"
                "   nnzPy = 0;\n"
                "   nnzMax = 0;\n"
                "   for (i = 0; i < " << m << "; i++) {\n"
                "      if (py[i] != 0.0) {\n"
                "         nnzPy++;\n"
                "         pyPos = (unsigned long*) realloc(pyPos, nnzPy * sizeof(unsigned long));\n"
                "         pyPos[nnzPy - 1] = i;\n"
                "         " << _name << "_" << FUNCTION_REVERSE_ONE_SPARSITY << "(i, &pos, &nnz);\n"
                "         if(nnz > nnzMax)\n"
                "            nnzMax = nnz;\n"
                "      }\n"
                "   }\n"
                "   for (j = 0; j < " << n << "; j++) {\n"
                "      px[j] = 0;\n"
                "   }\n"
                "\n"
                "   if (nnzPy == 0) {\n"
                "      free(pyPos);\n"
                "      return 0; //nothing to do\n"
                "   }\n"
                "\n"
                "   compressed = (" << _baseTypeName << "*) malloc(nnzMax * sizeof(" << _baseTypeName << "));\n"
                "\n"
                "   for (ei = 0; ei < nnzPy; ei++) {\n"
                "      i = pyPos[ei];\n"
                "      " << _name << "_" << FUNCTION_REVERSE_ONE_SPARSITY << "(i, &pos, &nnz);\n"
                "\n"
                "      in[0] = x;\n"
                "      in[1] = &py[i];\n"
                "      out[0] = compressed;\n"
                "      " << _name << "_" << FUNCTION_SPARSE_REVERSE_ONE << "(i, " << args << ");\n"
                "\n"
                "      for (ePos = 0; ePos < nnz; ePos++) {\n"
                "         px[pos[ePos]] += compressed[ePos];\n"
                "      }\n"
                "\n"
                "   }\n"
                "   free(compressed);\n"
                "   free(pyPos);\n"
                "   return 0;\n"
                "}\n";
        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseReverseTwoSources(std::map<std::string, std::string>& sources) throw (CGException) {
        const size_t m = _fun.Range();
        const size_t n = _fun.Domain();
        //const size_t k = 1;
        const size_t p = 2;

        determineHessianSparsity();

        /**
         * we might have to consider a slightly different order than the one
         * specified by the user according to the available elements in the sparsity
         */
        std::vector<size_t> evalRows, evalCols;
        determineSecondOrderElements4Eval(evalRows, evalCols);

        // elements[var]{vars}
        std::map<size_t, std::vector<size_t> > elements;
        for (size_t e = 0; e < evalCols.size(); e++) {
            elements[evalRows[e]].push_back(evalCols[e]);
        }

        if (!_loopTapes.empty()) {
            /**
             * with loops
             */
            prepareSparseReverseTwoWithLoops(sources, elements);
            return;
        }

        vector<CGBase> tx1v(n);

        std::set<size_t> customVarsInHess;
        if (_custom_hess.defined) {
            customVarsInHess.insert(_custom_hess.row.begin(), _custom_hess.row.end());
            customVarsInHess.insert(_custom_hess.col.begin(), _custom_hess.col.end());
        }

#if 0  // coloring must be changed
        /**
         * Coloring
         */
        printSparsityPattern(_hessSparsity.sparsity, "hessian");
        const CppAD::vector<Color> colors = colorByRow(customVarsInHess, _hessSparsity.sparsity);

        for (size_t c = 0; c < colors.size(); c++) {
            const Color& color = colors[c];

            _cache.str("");
            _cache << "model (reverse two, indep " << *color.rows.begin();
            std::set<size_t>::const_iterator iti = color.rows.begin();
            size_t count = 1;
            for (++iti; iti != color.rows.end() && count < 3; ++iti, count++) {
                _cache << " + " << *iti;
            }
            if (color.rows.size() > count) {
                _cache << *color.rows.begin() << " + {" << (color.rows.size() - count) << " others})";
            }
            _cache << ")";

            const std::string jobName = _cache.str();

            startingGraphCreation(jobName);

            CodeHandler<Base> handler;
            handler.setVerbose(_verbose);

            vector<CGBase> tx0(n);
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

            vector<CGBase> py(m); // (k+1)*m is not used because we are not interested in all values
            handler.makeVariables(py);
            if (_x.size() > 0) {
                for (size_t i = 0; i < m; i++) {
                    py[i].setValue(Base(1.0));
                }
            }

            _fun.Forward(0, tx0);

            for (std::set<size_t>::const_iterator itj = color.forbiddenRows.begin(); itj != color.forbiddenRows.end(); ++itj) {
                size_t j = *itj;
                tx1v[j] = tx1;
            }
            _fun.Forward(1, tx1v);

            // restore values
            for (std::set<size_t>::const_iterator itj = color.forbiddenRows.begin(); itj != color.forbiddenRows.end(); ++itj) {
                size_t j = *itj;
                tx1v[j] = Base(0);
            }

            vector<CGBase> px = _fun.Reverse(2, py);
            assert(px.size() == 2 * n);

            finishedGraphCreation();


            std::map<size_t, std::set<size_t> >::const_iterator itr2c;
            for (itr2c = color.row2Columns.begin(); itr2c != color.row2Columns.end(); ++itr2c) {
                size_t j = itr2c->first;
                const std::set<size_t>& cols = itr2c->second;

                vector<CGBase> pxCustom(cols.size());
                std::set<size_t>::const_iterator it2;
                size_t e = 0;
                for (it2 = cols.begin(); it2 != cols.end(); ++it2, e++) {
                    size_t jj = *it2;
                    pxCustom[e] = px[jj * p + 1]; // not interested in all values
                }

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
#endif

        /**
         * Generate one function for each independent variable
         */
        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t j = it->first;
            const std::vector<size_t>& cols = it->second;

            _cache.str("");
            _cache << "model (reverse two, indep " << j << ")";
            const std::string jobName = _cache.str();

            startingGraphCreation(jobName);

            CodeHandler<Base> handler;
            handler.setVerbose(_verbose);

            vector<CGBase> tx0(n);
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

            vector<CGBase> py(m); // (k+1)*m is not used because we are not interested in all values
            handler.makeVariables(py);
            if (_x.size() > 0) {
                for (size_t i = 0; i < m; i++) {
                    py[i].setValue(Base(1.0));
                }
            }

            _fun.Forward(0, tx0);

            tx1v[j] = tx1;
            _fun.Forward(1, tx1v);
            tx1v[j] = Base(0);
            vector<CGBase> px = _fun.Reverse(2, py);
            assert(px.size() == 2 * n);

            vector<CGBase> pxCustom;
            std::vector<size_t>::const_iterator it2;
            for (it2 = cols.begin(); it2 != cols.end(); ++it2) {
                size_t jj = *it2;
                pxCustom.push_back(px[jj * p + 1]); // not interested in all values
            }

            finishedGraphCreation();

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

        _cache.str("");

        generateGlobalDirectionalFunctionSource(FUNCTION_SPARSE_REVERSE_TWO,
                                                "indep",
                                                FUNCTION_REVERSE_TWO_SPARSITY,
                                                elements,
                                                sources);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateReverseTwoSources(std::map<std::string, std::string>& sources) {
        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        _cache.str("");
        _cache << _name << "_" << FUNCTION_REVERSE_TWO;
        std::string model_function(_cache.str());
        _cache.str("");

        CLanguage<Base> langC(_baseTypeName);
        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();
        std::string args = langC.generateDefaultFunctionArguments();

        _cache << "#include <stdlib.h>\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n"
                "\n"
                "void " << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO << "(unsigned long pos, " << argsDcl << ");\n"
                "void " << _name << "_" << FUNCTION_REVERSE_TWO_SPARSITY << "(unsigned long pos, unsigned long const** elements, unsigned long* nnz);\n"
                "\n"
                "int " << model_function << "("
                << _baseTypeName << " const tx[], "
                << _baseTypeName << " const ty[], "
                << _baseTypeName << " px[], "
                << _baseTypeName << " const py[], "
                << langC.generateArgumentAtomicDcl() << ") {\n"
                "    unsigned long ej, ePos, i, j, nnz, nnzMax;\n"
                "    unsigned long const* pos;\n"
                "    unsigned long* txPos;\n"
                "    unsigned long nnzTx;\n"
                "    " << _baseTypeName << " const * in[3];\n"
                "    " << _baseTypeName << "* out[1];\n"
                "    " << _baseTypeName << " x[" << n << "];\n"
                "    " << _baseTypeName << " w[" << m << "];\n"
                "    " << _baseTypeName << "* compressed;\n"
                "    int nonZeroW;\n"
                "\n"
                "    nonZeroW = 0;\n"
                "    for (i = 0; i < " << m << "; i++) {\n"
                "        if (py[i * 2] != 0.0) {\n"
                "            return 1; // error\n"
                "        }\n"
                "        w[i] = py[i * 2 + 1];\n"
                "        if(w[i] != 0.0) nonZeroW++;\n"
                "    }\n"
                "\n"
                "    for (j = 0; j < " << n << "; j++) {\n"
                "       px[j * 2] = 0;\n"
                "    }\n"
                "\n"
                "    if(nonZeroW == 0)\n"
                "        return 0; //nothing to do\n"
                "\n"
                "   txPos = 0;\n"
                "   nnzTx = 0;\n"
                "   nnzMax = 0;\n"
                "   for (j = 0; j < " << n << "; j++) {\n"
                "      if (tx[j * 2 + 1] != 0.0) {\n"
                "         nnzTx++;\n"
                "         txPos = (unsigned long*) realloc(txPos, nnzTx * sizeof(unsigned long));\n"
                "         txPos[nnzTx - 1] = j;\n"
                "         " << _name << "_" << FUNCTION_REVERSE_TWO_SPARSITY << "(j, &pos, &nnz);\n"
                "         if(nnz > nnzMax)\n"
                "            nnzMax = nnz;\n"
                "      }\n"
                "   }\n"
                "\n"
                "   if (nnzTx == 0) {\n"
                "      free(txPos);\n"
                "      return 0; //nothing to do\n"
                "   }\n"
                "\n"
                "    for (j = 0; j < " << n << "; j++)\n"
                "        x[j] = tx[j * 2];\n"
                "\n"
                "   compressed = (" << _baseTypeName << "*) malloc(nnzMax * sizeof(" << _baseTypeName << "));\n"
                "\n"
                "   for (ej = 0; ej < nnzTx; ej++) {\n"
                "      j = txPos[ej];\n"
                "      " << _name << "_" << FUNCTION_REVERSE_TWO_SPARSITY << "(j, &pos, &nnz);\n"
                "\n"
                "      in[0] = x;\n"
                "      in[1] = &tx[j * 2 + 1];\n"
                "      in[2] = w;\n"
                "      out[0] = compressed;\n";
        if (!_loopTapes.empty()) {
            _cache << "      for(ePos = 0; ePos < nnz; ePos++)\n"
                    "         compressed[ePos] = 0;\n"
                    "\n";
        }
        _cache << "      " << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO << "(j, " << args << ");\n"
                "\n"
                "      for (ePos = 0; ePos < nnz; ePos++) {\n"
                "         px[pos[ePos] * 2] += compressed[ePos];\n"
                "      }\n"
                "\n"
                "   }\n"
                "   free(compressed);\n"
                "   free(txPos);\n"
                "   return 0;\n"
                "};\n";

        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateGlobalDirectionalFunctionSource(const std::string& function,
                                                                                const std::string& suffix,
                                                                                const std::string& function_sparsity,
                                                                                const std::map<size_t, std::vector<size_t> >& elements,
                                                                                std::map<std::string, std::string>& sources) {
        /**
         * The function that matches each equation to a directional derivative function
         */
        CLanguage<Base> langC(_baseTypeName);
        std::string argsDcl = langC.generateDefaultFunctionArgumentsDcl();
        std::string args = langC.generateDefaultFunctionArguments();

        _cache.str("");
        _cache << _name << "_" << function;
        std::string model_function = _cache.str();
        _cache.str("");

        _cache << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n\n";
        generateFunctionDeclarationSource(_cache, model_function, suffix, elements, argsDcl);
        _cache << "\n";
        _cache << "int " << model_function << "("
                "unsigned long pos, " << argsDcl << ") {\n"
                "   switch(pos) {\n";
        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            // the size of each sparsity row
            _cache << "      case " << it->first << ":\n"
                    "         " << model_function << "_" << suffix << it->first << "(" << args << ");\n"
                    "         return 0; // done\n";
        }
        _cache << "      default:\n"
                "         return 1; // error\n"
                "   };\n";

        _cache << "}\n";
        sources[model_function + ".c"] = _cache.str();
        _cache.str("");

        /**
         * Sparsity
         */
        generateSparsity1DSource2(_name + "_" + function_sparsity, elements);
        sources[_name + "_" + function_sparsity + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    template<class T>
    void CLangCompileModelHelper<Base>::generateFunctionDeclarationSource(std::ostringstream& cache,
                                                                          const std::string& model_function,
                                                                          const std::string& suffix,
                                                                          const std::map<size_t, T>& elements,
                                                                          const std::string& argsDcl) {
        typename std::map<size_t, T>::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t pos = it->first;
            cache << "void " << model_function << "_" << suffix << pos << "(" << argsDcl << ");\n";
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparsity1DSource(const std::string& function,
                                                                 const std::vector<size_t>& sparsity) {
        _cache << "void " << function << "("
                "unsigned long const** sparsity,"
                " unsigned long* nnz) {\n";

        // the size of each sparsity row
        _cache << "   static unsigned long const nonzeros[" << sparsity.size() << "] = {";
        if (!sparsity.empty()) {
            _cache << sparsity[0];
            for (size_t i = 1; i < sparsity.size(); i++) {
                _cache << "," << sparsity[i];
            }
        }
        _cache << "};\n";

        _cache << "   *sparsity = nonzeros;\n"
                "   *nnz = " << sparsity.size() << ";\n"
                "}\n";
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparsity2DSource(const std::string& function,
                                                                 const LocalSparsityInfo& sparsity) {
        const std::vector<size_t>& rows = sparsity.rows;
        const std::vector<size_t>& cols = sparsity.cols;

        assert(rows.size() == cols.size());

        _cache << "void " << function << "("
                "unsigned long const** row,"
                " unsigned long const** col,"
                " unsigned long* nnz) {\n";

        // the size of each sparsity row
        _cache << "static unsigned long const rows[" << rows.size() << "] = {";
        if (!rows.empty()) {
            _cache << rows[0];
            for (size_t i = 1; i < rows.size(); i++) {
                _cache << "," << rows[i];
            }
        }
        _cache << "};\n";

        _cache << "static unsigned long const cols[" << cols.size() << "] = {";
        if (!cols.empty()) {
            _cache << cols[0];
            for (size_t i = 1; i < cols.size(); i++) {
                _cache << "," << cols[i];
            }
        }
        _cache << "};\n";

        _cache << "*row = rows;\n"
                "*col = cols;\n"
                "*nnz = " << rows.size() << ";\n"
                "}\n";
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparsity2DSource2(const std::string& function,
                                                                  const std::vector<LocalSparsityInfo>& sparsities) {
        _cache << "void " << function << "("
                "unsigned long i,"
                "unsigned long const** row,"
                " unsigned long const** col,"
                " unsigned long* nnz) {\n";

        for (size_t i = 0; i < sparsities.size(); i++) {
            const std::vector<size_t>& rows = sparsities[i].rows;
            const std::vector<size_t>& cols = sparsities[i].cols;
            assert(rows.size() == cols.size());
            if (!rows.empty()) {
                _cache << "   static unsigned long const rows" << i << "[" << rows.size() << "] = {";
                _cache << rows[0];
                for (size_t j = 1; j < rows.size(); j++) {
                    _cache << "," << rows[j];
                }
                _cache << "};\n";

                _cache << "   static unsigned long const cols" << i << "[" << cols.size() << "] = {";
                _cache << cols[0];
                for (size_t i = 1; i < cols.size(); i++) {
                    _cache << "," << cols[i];
                }
                _cache << "};\n";
            }
        }

        _cache << "   switch(i) {\n";
        for (size_t i = 0; i < sparsities.size(); i++) {
            // the size of each sparsity
            if (!sparsities[i].rows.empty()) {
                _cache << "   case " << i << ":\n"
                        "      *row = rows" << i << ";\n"
                        "      *col = cols" << i << ";\n"
                        "      *nnz = " << sparsities[i].rows.size() << ";\n"
                        "      break;\n";
            }
        }

        _cache << "   default:\n"
                "      *row = 0;\n"
                "      *col = 0;\n"
                "      *nnz = 0;\n"
                "   break;\n"
                "   };\n"
                "}\n";
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparsity1DSource2(const std::string& function,
                                                                  const std::map<size_t, std::vector<size_t> >& elements) {

        _cache << "void " << function << "("
                "unsigned long pos,"
                " unsigned long const** elements,"
                " unsigned long* nnz) {\n";

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            // the size of each sparsity row
            const std::vector<size_t>& els = it->second;
            _cache << "   static unsigned long const elements" << it->first << "[" << els.size() << "] = {";
            if (!els.empty()) {
                _cache << els[0];
                for (size_t i = 1; i < els.size(); i++) {
                    _cache << "," << els[i];
                }
            }
            _cache << "};\n";
        }

        _cache << "   switch(pos) {\n";
        for (it = elements.begin(); it != elements.end(); ++it) {
            // the size of each sparsity row
            _cache << "   case " << it->first << ":\n"
                    "      *elements = elements" << it->first << ";\n"
                    "      *nnz = " << it->second.size() << ";\n"
                    "      break;\n";
        }
        _cache << "   default:\n"
                "      *elements = 0;\n"
                "      *nnz = 0;\n"
                "   break;\n"
                "   };\n"
                "}\n";
    }

    template<class Base>
    void inline CLangCompileModelHelper<Base>::startingGraphCreation(const std::string& jobName) {
        if (_verbose) {
            if (!_jobNames.empty()) {
                if (!_nestedJobs.back())
                    std::cout << "\n";
                std::string ident(3 * _jobNames.size(), ' ');
                std::cout << ident;
            }
            std::cout << "generating operation graph for '" << jobName << "' ... ";
            std::cout.flush();

            std::fill(_nestedJobs.begin(), _nestedJobs.end(), true);

            _nestedJobs.push_back(false);
            _jobNames.push_back(jobName);
            _beginTimes.push_back(system::currentTime());
        }
    }

    template<class Base>
    void inline CLangCompileModelHelper<Base>::finishedGraphCreation() {
        if (_verbose) {
            assert(_nestedJobs.size() > 0);

            double beginTime = _beginTimes.back();
            double endTime = system::currentTime();
            double elapsed = endTime - beginTime;
            if (_nestedJobs.back()) {
                std::string ident(3 * (_jobNames.size() - 1), ' ');
                std::cout << ident;
                std::cout << "generated operation graph for '" << _jobNames.back() << "' ... ";
            }

            std::cout << "done [" << std::fixed << std::setprecision(3) << elapsed << "]" << std::endl;

            _nestedJobs.pop_back();
            _jobNames.pop_back();
            _beginTimes.pop_back();
        }
    }

    template<class Base>
    inline std::map<size_t, std::vector<std::set<size_t> > > CLangCompileModelHelper<Base>::determineOrderByCol(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                                                const LocalSparsityInfo& sparsity) {
        return determineOrderByCol(elements, sparsity.rows, sparsity.cols);
    }

    template<class Base>
    inline std::map<size_t, std::vector<std::set<size_t> > > CLangCompileModelHelper<Base>::determineOrderByCol(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                                                const std::vector<size_t>& userRows,
                                                                                                                const std::vector<size_t>& userCols) {
        std::map<size_t, std::vector<std::set<size_t> > > userLocation;

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t col = it->first;
            const std::vector<size_t>& els = it->second;
            std::vector<std::set<size_t> >& userLocationCol = userLocation[col];
            userLocationCol.resize(els.size());

            for (size_t er = 0; er < els.size(); er++) {
                size_t row = els[er];
                for (size_t e = 0; e < userRows.size(); e++) {
                    if (userRows[e] == row && userCols[e] == col) {
                        userLocationCol[er].insert(e);
                        break;
                    }
                }
            }
        }

        return userLocation;
    }

    template<class Base>
    inline std::map<size_t, std::vector<std::set<size_t> > > CLangCompileModelHelper<Base>::determineOrderByRow(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                                                const LocalSparsityInfo& sparsity) {
        return determineOrderByRow(elements, sparsity.rows, sparsity.cols);
    }

    template<class Base>
    inline std::map<size_t, std::vector<std::set<size_t> > > CLangCompileModelHelper<Base>::determineOrderByRow(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                                                const std::vector<size_t>& userRows,
                                                                                                                const std::vector<size_t>& userCols) {
        std::map<size_t, std::vector<std::set<size_t> > > userLocation;

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t row = it->first;
            const std::vector<size_t>& els = it->second;
            std::vector<std::set<size_t> >& userLocationRow = userLocation[row];
            userLocationRow.resize(els.size());

            for (size_t ec = 0; ec < els.size(); ec++) {
                size_t col = els[ec];
                for (size_t e = 0; e < userRows.size(); e++) {
                    if (userCols[e] == col && userRows[e] == row) {
                        userLocationRow[ec].insert(e);
                        break;
                    }
                }
            }
        }

        return userLocation;
    }

    /**
     * 
     * Specializations
     */
    template<>
    inline std::string CLangCompileModelHelper<double>::baseTypeName() {
        return "double";
    }

    template<>
    inline std::string CLangCompileModelHelper<float>::baseTypeName() {
        return "float";
    }
}

#endif
