#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_IMPL_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_IMPL_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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
    const std::string CLangCompileModelHelper<Base>::CONST = "const";

    template<class Base>
    VariableNameGenerator<Base>* CLangCompileModelHelper<Base>::createVariableNameGenerator(const std::string& depName,
                                                                                            const std::string& indepName,
                                                                                            const std::string& tmpName) {
        return new CLangDefaultVariableNameGenerator<Base > (depName, indepName, tmpName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::compileSources(CLangCompiler<Base>& compiler,
                                                       bool posIndepCode) {
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

        if (_sparseJacobian) {
            generateSparseJacobianSource(sources);
        }

        if (_sparseHessian) {
            generateSparseHessianSource(sources);
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

        if (_sparseJacobian || _forwardOne || _reverseOne) {
            generateJacobianSparsitySource(sources);
        }

        if (_sparseHessian || _reverseTwo) {
            generateHessianSparsitySource(sources);
        }

        generateInfoSource(sources);

        compiler.compileSources(sources, posIndepCode, true);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateInfoSource(std::map<std::string, std::string>& sources) {
        const char* localBaseName = typeid (Base).name();

        std::string funcName = _name + "_" + FUNCTION_INFO;

        std::auto_ptr<VariableNameGenerator< Base > > nameGen(createVariableNameGenerator("dep", "ind", "var"));

        _cache.str("");
        _cache << "void " << funcName << "(const char** baseName, unsigned long int* m, unsigned long int* n, unsigned int* indCount, unsigned int* depCount) {\n"
                "   *baseName = \"" << _baseTypeName << "  " << localBaseName << "\";\n"
                "   *m = " << _fun->Range() << ";\n"
                "   *n = " << _fun->Domain() << ";\n"
                "   *depCount = " << nameGen->getDependent().size() << "; // number of dependent array variables\n"
                "   *indCount = " << nameGen->getIndependent().size() << "; // number of independent array variables\n"
                "}\n\n";

        sources[funcName + ".c"] = _cache.str();
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateZeroSource(std::map<std::string, std::string>& sources) {
        const std::string jobName = "model (zero-order forward)";

        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        std::vector<CGBase> indVars(_fun->Domain());
        handler.makeVariables(indVars);

        std::vector<CGBase> dep = _fun->Forward(0, indVars);

        finishedGraphCreation();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_FORWAD_ZERO);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("dep", "ind", "var"));

        handler.generateCode(code, langC, dep, *nameGen, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateJacobianSource(std::map<std::string, std::string>& sources) {
        const std::string jobName = "Jacobian";

        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        std::vector<CGBase> indVars(_fun->Domain());
        handler.makeVariables(indVars);

        std::vector<CGBase> jac = _fun->Jacobian(indVars);

        finishedGraphCreation();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_JACOBIAN);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("jac", "ind", "var"));

        handler.generateCode(code, langC, jac, *nameGen, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateHessianSource(std::map<std::string, std::string>& sources) {
        const std::string jobName = "Hessian";

        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        //size_t m = _fun->Range();
        size_t n = _fun->Domain();


        // independent variables
        std::vector<CGBase> indVars(n);
        handler.makeVariables(indVars);
        // multipliers
        std::vector<CGBase> w(_fun->Range());
        handler.makeVariables(w);

        std::vector<CGBase> hess = _fun->Hessian(indVars, w);

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
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("hess", "ind", "var"));
        CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), n);

        handler.generateCode(code, langC, hess, nameGenHess, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseJacobianSource(std::map<std::string, std::string>& sources) {
        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        /**
         * Determine the sparsity pattern
         */
        determineJacobianSparsity();

        /**
         * Estimate the work load of forward vs reverse mode
         */
        size_t workForward;
        size_t workReverse;
        if (_custom_jac.defined) {
            std::set<size_t> rows, cols;
            rows.insert(_jacSparsity.rows.begin(), _jacSparsity.rows.end());
            workReverse = rows.size();
            cols.insert(_jacSparsity.cols.begin(), _jacSparsity.cols.end());
            workForward = cols.size();
        } else {
            workReverse = m;
            workForward = n;
        }

        /**
         * call the appropriate method for source code generation
         */
        if (_forwardOne && workForward <= workReverse) {
            generateSparseJacobianForRevSource(sources, true);
        } else if (_reverseOne && workForward >= workReverse) {
            generateSparseJacobianForRevSource(sources, false);
        } else {
            generateSparseJacobianSource(sources, workForward <= workReverse);
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseJacobianSource(std::map<std::string, std::string>& sources,
                                                                     bool forward) {
        const std::string jobName = "sparse Jacobian";

        //size_t m = _fun->Range();
        size_t n = _fun->Domain();

        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        std::vector<CGBase> indVars(n);
        handler.makeVariables(indVars);

        std::vector<CGBase> jac(_jacSparsity.rows.size());
        CppAD::sparse_jacobian_work work;
        if (forward) {
            _fun->SparseJacobianForward(indVars, _jacSparsity.sparsity, _jacSparsity.rows, _jacSparsity.cols, jac, work);
        } else {
            _fun->SparseJacobianReverse(indVars, _jacSparsity.sparsity, _jacSparsity.rows, _jacSparsity.cols, jac, work);
        }

        finishedGraphCreation();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_SPARSE_JACOBIAN);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("jac", "ind", "var"));

        handler.generateCode(code, langC, jac, *nameGen, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseJacobianForRevSource(std::map<std::string, std::string>& sources,
                                                                           bool forward) {
        //size_t m = _fun->Range();
        //size_t n = _fun->Domain();
        using namespace std;
        _cache.str("");
        _cache << _name << "_" << FUNCTION_SPARSE_JACOBIAN;
        string model_function(_cache.str());

        map<size_t, std::vector<size_t> >::const_iterator it;

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

        size_t maxCompressedSize = 0;
        for (it = elements.begin(); it != elements.end(); ++it) {
            if (it->second.size() > maxCompressedSize)
                maxCompressedSize = it->second.size();
        }

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                "\n";
        generateFunctionDeclarationSource(_cache, functionRevFor, revForSuffix, elements);
        _cache << "\n"
                "void " << model_function << "(" << _baseTypeName << " const *const * in, " << _baseTypeName << " *const * out) {\n"
                "   " << _baseTypeName << " const * inLocal[2];\n"
                "   " << _baseTypeName << " inLocal1 = 1;\n"
                "   " << _baseTypeName << " * outLocal[1];\n"
                "   " << _baseTypeName << " compressed[" << maxCompressedSize << "];\n"
                "   " << _baseTypeName << " * jac = out[0];\n"
                "\n"
                "   inLocal[0] = in[0];\n"
                "   inLocal[1] = &inLocal1;\n"
                "   outLocal[0] = compressed;";
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t index = it->first;
            const std::vector<size_t>& els = it->second;
            const std::vector<set<size_t> >& location = userJacElLocation.at(index);
            assert(els.size() == location.size());

            _cache << "\n"
                    "   " << functionRevFor << "_" << revForSuffix << index << "(inLocal, outLocal);\n";
            for (size_t e = 0; e < els.size(); e++) {
                _cache << "   ";
                set<size_t>::const_iterator itl;
                for(itl = location[e].begin(); itl != location[e].end(); ++itl) {
                    _cache << "jac[" << (*itl) << "] = ";
                }
                _cache << "compressed[" << e << "];\n";
            }
        }

        _cache << "}\n";
        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::determineJacobianSparsity() {
        if (!_jacSparsity.sparsity.empty()) {
            return;
        }
        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        /**
         * Determine the sparsity pattern
         */
        _jacSparsity.sparsity = jacobianSparsity < std::vector<bool>, CGBase > (*_fun);

        if (!_custom_jac.defined) {
            generateSparsityIndexes(_jacSparsity.sparsity, m, n, _jacSparsity.rows, _jacSparsity.cols);

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
    void CLangCompileModelHelper<Base>::generateSparseHessianSource(std::map<std::string, std::string>& sources) {
        const std::string jobName = "sparse Hessian";
        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        /**
         * Determine the sparsity pattern p for Hessian of w^T F
         */
        determineHessianSparsity();


        // make use of the symmetry of the Hessian in order to reduce operations
        std::map<size_t, std::map<size_t, size_t> > locations;
        for (size_t i = 0; i < _hessSparsity.rows.size(); i++) {
            locations[_hessSparsity.rows[i]][_hessSparsity.cols[i]] = i;
        }

        std::vector<size_t> upperHessRows, upperHessCols, upperHessOrder;
        upperHessRows.reserve(_hessSparsity.rows.size() / 2);
        upperHessCols.reserve(upperHessRows.size());
        upperHessOrder.reserve(upperHessRows.size());

        std::map<size_t, size_t> duplicates; // the elements determined using symmetry
        std::map<size_t, std::map<size_t, size_t> >::const_iterator ii;
        std::map<size_t, size_t>::const_iterator jj;
        for (size_t i = 0; i < _hessSparsity.rows.size(); i++) {
            bool add = true;
            if (_hessSparsity.rows[i] > _hessSparsity.cols[i]) {
                ii = locations.find(_hessSparsity.cols[i]);
                if (ii != locations.end()) {
                    jj = ii->second.find(_hessSparsity.rows[i]);
                    if (jj != ii->second.end()) {
                        size_t k = jj->second;
                        duplicates[i] = k;
                        add = false; // symmetric value being determined
                    }
                }
            }

            if (add) {
                upperHessRows.push_back(_hessSparsity.rows[i]);
                upperHessCols.push_back(_hessSparsity.cols[i]);
                upperHessOrder.push_back(i);
            }
        }

        /**
         * 
         */
        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        // independent variables
        std::vector<CGBase> indVars(n);
        handler.makeVariables(indVars);
        // multipliers
        std::vector<CGBase> w(m);
        handler.makeVariables(w);

        CppAD::sparse_hessian_work work;
        std::vector<CGBase> upperHess(upperHessRows.size());
        _fun->SparseHessian(indVars, w, _hessSparsity.sparsity, upperHessRows, upperHessCols, upperHess, work);

        std::vector<CGBase> hess(_hessSparsity.rows.size());
        for (size_t i = 0; i < upperHessOrder.size(); i++) {
            hess[upperHessOrder[i]] = upperHess[i];
        }

        // make use of the symmetry of the Hessian in order to reduce operations
        std::map<size_t, size_t>::const_iterator it2;
        for (it2 = duplicates.begin(); it2 != duplicates.end(); ++it2) {
            hess[it2->first] = hess[it2->second];
        }

        finishedGraphCreation();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_SPARSE_HESSIAN);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("hess", "ind", "var"));
        CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), n);

        handler.generateCode(code, langC, hess, nameGenHess, jobName);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::determineHessianSparsity() {
        if (!_hessSparsity.sparsity.empty()) {
            return;
        }

        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        _hessSparsity.sparsity = hessianSparsity < std::vector<bool>, CGBase > (*_fun);

        if (!_custom_hess.defined) {
            generateSparsityIndexes(_hessSparsity.sparsity, n, n,
                                    _hessSparsity.rows, _hessSparsity.cols);

        } else {
            _hessSparsity.rows = _custom_hess.row;
            _hessSparsity.cols = _custom_hess.col;
        }

        /**
         * For each individual equation
         */
        _hessSparsities.resize(m);
        for (size_t i = 0; i < m; i++) {
            _hessSparsities[i].sparsity = hessianSparsity < std::vector<bool>, CGBase > (*_fun, i);

            if (!_custom_hess.defined) {
                generateSparsityIndexes(_hessSparsities[i].sparsity, n, n,
                                        _hessSparsities[i].rows, _hessSparsities[i].cols);

            } else {
                for (size_t e = 0; e < _custom_hess.row.size(); e++) {
                    size_t i1 = _custom_hess.row[e];
                    size_t i2 = _custom_hess.col[e];
                    if (_hessSparsities[i].sparsity[i1 * n + i2]) {
                        _hessSparsities[i].rows.push_back(i1);
                        _hessSparsities[i].cols.push_back(i2);
                    }
                }
            }
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateHessianSparsitySource(std::map<std::string, std::string>& sources) {
        determineHessianSparsity();

        generateSparsity2DSource(_name + "_" + FUNCTION_HESSIAN_SPARSITY, _hessSparsity);
        sources[_name + "_" + FUNCTION_HESSIAN_SPARSITY + ".c"] = _cache.str();
        _cache.str("");

        generateSparsity2DSource2(_name + "_" + FUNCTION_HESSIAN_SPARSITY2, _hessSparsities);
        sources[_name + "_" + FUNCTION_HESSIAN_SPARSITY2 + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseForwardOneSources(std::map<std::string, std::string>& sources) {
        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        determineJacobianSparsity();

        // elements[var]{equations}
        std::map<size_t, std::vector<size_t> > elements;
        for (size_t e = 0; e < _jacSparsity.rows.size(); e++) {
            elements[_jacSparsity.cols[e]].push_back(_jacSparsity.rows[e]);
        }

        /**
         * Generate one function for each dependent variable
         */
        std::vector<CGBase> dxv(n, Base(0));

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

            std::vector<CGBase> indVars(n);
            handler.makeVariables(indVars);

            CGBase dx;
            handler.makeVariable(dx);

            // TODO: consider caching the zero order coefficients somehow between calls
            _fun->Forward(0, indVars);
            dxv[j] = dx;
            std::vector<CGBase> dy = _fun->Forward(1, dxv);
            dxv[j] = Base(0);
            assert(dy.size() == m);

            std::vector<CGBase> dyCustom;
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
            std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("dy", "ind", "var"));
            CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), "dx", n);

            handler.generateCode(code, langC, dyCustom, nameGenHess, jobName);
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

        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        _cache.str("");
        _cache << _name << "_" << FUNCTION_FORWARD_ONE;
        std::string model_function(_cache.str());

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                "\n"
                "void " << _name << "_" << FUNCTION_SPARSE_FORWARD_ONE << "(unsigned long int pos," << _baseTypeName << " const *const * in, " << _baseTypeName << " *const * out);\n"
                "void " << _name << "_" << FUNCTION_FORWARD_ONE_SPARSITY << "(unsigned long int pos, unsigned long int const** elements, unsigned long int* nnz);\n"
                "\n"
                "int " << model_function << "(" << _baseTypeName << " const tx[], " << _baseTypeName << " ty[]) {\n"
                "   unsigned long int e, i, j, jj, nnz;\n"
                "   unsigned long int const* pos;\n"
                "   " << _baseTypeName << " const * in[2];\n"
                "   " << _baseTypeName << "* out[1];\n"
                "   " << _baseTypeName << " x[" << n << "];\n"
                "   " << _baseTypeName << "* compressed;\n"
                "   int found;\n"
                "\n"
                "   found = 0;\n"
                "   for (jj = 0; jj < " << n << "; jj++) {\n"
                "      if (tx[jj * 2 + 1] != 0.0) {\n"
                "         if (found) {\n"
                "            return 1; // error \n"
                "         }\n"
                "         j = jj;\n"
                "         found = 1;\n"
                "      }\n"
                "   }\n"
                "   for (i = 0; i < " << m << "; i++) {\n"
                "      ty[i * 2 + 1] = 0;\n"
                "   }\n"
                "   if (!found) {\n"
                "      return 0; //nothing to do\n"
                "   }\n"
                "\n"
                "   " << _name << "_" << FUNCTION_FORWARD_ONE_SPARSITY << "(j, &pos, &nnz);\n"
                "\n"
                "   for (jj = 0; jj < " << n << "; jj++)\n"
                "      x[jj] = tx[jj * 2];\n"
                "\n"
                "   compressed = (" << _baseTypeName << "*) malloc(nnz * sizeof(" << _baseTypeName << "));\n"
                "   in[0] = x;\n"
                "   in[1] = &tx[j * 2 + 1];\n"
                "   out[0] = compressed;\n"
                "   " << _name << "_" << FUNCTION_SPARSE_FORWARD_ONE << "(j, in, out);\n"
                "\n"
                "   for (e = 0; e < nnz; e++) {\n"
                "      ty[pos[e] * 2 + 1] = compressed[e];\n"
                "   }\n"
                "\n"
                "   free(compressed);\n"
                "   return 0;\n"
                "}\n";
        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseReverseOneSources(std::map<std::string, std::string>& sources) {
        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        determineJacobianSparsity();

        // elements[equation]{vars}
        std::map<size_t, std::vector<size_t> > elements;
        for (size_t e = 0; e < _jacSparsity.rows.size(); e++) {
            elements[_jacSparsity.rows[e]].push_back(_jacSparsity.cols[e]);
        }

        std::vector<CGBase> w(m, Base(0));

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

            std::vector<CGBase> indVars(_fun->Domain());
            handler.makeVariables(indVars);

            CGBase py;
            handler.makeVariable(py);

            // TODO: consider caching the zero order coefficients somehow between calls
            _fun->Forward(0, indVars);

            w[i] = py;
            std::vector<CGBase> dw = _fun->Reverse(1, w);
            assert(dw.size() == n);
            w[i] = Base(0);

            std::vector<CGBase> dwCustom;
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
            std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("dw", "ind", "var"));
            CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), "py", n);

            handler.generateCode(code, langC, dwCustom, nameGenHess, jobName);
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
        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        _cache.str("");
        _cache << _name << "_" << FUNCTION_REVERSE_ONE;
        std::string model_function(_cache.str());

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                "\n"
                "void " << _name << "_" << FUNCTION_SPARSE_REVERSE_ONE << "(unsigned long int pos," << _baseTypeName << " const *const * in, " << _baseTypeName << " *const * out);\n"
                "void " << _name << "_" << FUNCTION_REVERSE_ONE_SPARSITY << "(unsigned long int pos, unsigned long int const** elements, unsigned long int* nnz);\n"
                "\n"
                "int " << model_function << "("
                << _baseTypeName << " const x[], "
                << _baseTypeName << " const ty[],"
                << _baseTypeName << "  px[], "
                << _baseTypeName << " const py[]) {\n"
                "   unsigned long int e, i, ii, j, nnz;\n"
                "   unsigned long int const* pos;\n"
                "   " << _baseTypeName << " const * in[2];\n"
                "   " << _baseTypeName << "* out[1];\n"
                "   " << _baseTypeName << "* compressed;\n"
                "   int found;\n"
                "\n"
                "   found = 0;\n"
                "   for (ii = 0; ii < " << m << "; ii++) {\n"
                "      if (py[ii] != 0.0) {\n"
                "         if (found) {\n"
                "            return 1; // error \n"
                "         }\n"
                "         i = ii;\n"
                "         found = 1;\n"
                "      }\n"
                "   }\n"
                "   for (j = 0; j < " << n << "; j++) {\n"
                "      px[j] = 0;\n"
                "   }\n"
                "   if (!found) {\n"
                "      return 0; //nothing to do\n"
                "   }\n"
                "\n"
                "   " << _name << "_" << FUNCTION_REVERSE_ONE_SPARSITY << "(i, &pos, &nnz);\n"
                "\n"
                "   compressed = (" << _baseTypeName << "*) malloc(nnz * sizeof(" << _baseTypeName << "));\n"
                "   in[0] = x;\n"
                "   in[1] = &py[i];\n"
                "   out[0] = compressed;\n"
                "   " << _name << "_" << FUNCTION_SPARSE_REVERSE_ONE << "(i, in, out);\n"
                "\n"
                "   for (e = 0; e < nnz; e++) {\n"
                "      px[pos[e]] = compressed[e];\n"
                "   }\n"
                "\n"
                "   free(compressed);\n"
                "   return 0;\n"
                "}\n";
        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseReverseTwoSources(std::map<std::string, std::string>& sources) {
        const size_t m = _fun->Range();
        const size_t n = _fun->Domain();
        const size_t k = 1;
        const size_t k1 = k + 1;
        const size_t p = 2;

        determineHessianSparsity();

        // elements[vars]{equation}
        std::map<size_t, std::vector<size_t> > elements;
        for (size_t e = 0; e < _hessSparsity.cols.size(); e++) {
            elements[_hessSparsity.cols[e]].push_back(_hessSparsity.rows[e]);
        }

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

            CodeHandler<Base> handler;
            handler.setVerbose(_verbose);

            std::vector<CGBase> tx0(n);
            handler.makeVariables(tx0);

            CGBase tx1;
            handler.makeVariable(tx1);

            std::vector<CGBase> w(k1 * m);
            handler.makeVariables(w);

            // TODO: consider caching the zero order coefficients somehow between calls
            std::vector<CGBase> y = _fun->Forward(0, tx0);

            tx1v[j] = tx1;
            std::vector<CGBase> y_p = _fun->Forward(1, tx1v);
            tx1v[j] = Base(0);
            std::vector<CGBase> ddw = _fun->Reverse(2, w);
            assert(ddw.size() == 2 * n);

            std::vector<CGBase> ddwCustom;
            std::vector<size_t>::const_iterator it2;
            for (it2 = cols.begin(); it2 != cols.end(); ++it2) {
                size_t jj = *it2;
                ddwCustom.push_back(ddw[jj * p]);
            }

            finishedGraphCreation();

            CLanguage<Base> langC(_baseTypeName);
            langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
            _cache.str("");
            _cache << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO << "_indep" << j;
            langC.setGenerateFunction(_cache.str());

            std::ostringstream code;
            std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("px", "ind", "var"));
            CLangDefaultReverse2VarNameGenerator<Base> nameGenRev2(nameGen.get(), n, 1);

            handler.generateCode(code, langC, ddwCustom, nameGenRev2, jobName);
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
        //size_t m = _fun->Range();
        size_t n = _fun->Domain();

        _cache.str("");
        _cache << _name << "_" << FUNCTION_REVERSE_TWO;
        std::string model_function(_cache.str());

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                "\n"
                "void " << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO << "(unsigned long int pos, " << _baseTypeName << " const *const * in, " << _baseTypeName << " * const * out);\n"
                "void " << _name << "_" << FUNCTION_REVERSE_TWO_SPARSITY << "(unsigned long int pos, unsigned long int const** elements, unsigned long int* nnz);\n"
                "\n"
                "int " << model_function << "(" << _baseTypeName << " const tx[], " << _baseTypeName << " const ty[], " << _baseTypeName << " px[], " << _baseTypeName << " const py[]) {\n"
                "    unsigned long int e, i, j, jj, nnz;\n"
                "    unsigned long int const* pos;\n"
                "    " << _baseTypeName << " const * in[3];\n"
                "    " << _baseTypeName << " * out[1];\n"
                "    " << _baseTypeName << " x[" << n << "];\n"
                "    " << _baseTypeName << "* compressed;\n"
                "    int found;\n"
                "\n"
                "    found = 0;\n"
                "    for (jj = 0; jj < " << n << "; jj++) {\n"
                "        if (tx[jj * 2 + 1] != 0.0) {\n"
                "            if (found) {\n"
                "                return 1; // error \n"
                "            }\n"
                "            j = jj;\n"
                "            found = 1;\n"
                "        }\n"
                "    }\n"
                "    for (jj = 0; jj < " << n << "; jj++) {\n"
                "        px[jj * 2] = 0;\n"
                "    }\n"
                "    if (!found) {\n"
                "        return 0; //nothing to do\n"
                "    }\n"
                "\n"
                "    for (jj = 0; jj < " << n << "; jj++)\n"
                "        x[jj] = tx[jj * 2];\n"
                "\n"
                "    " << _name << "_" << FUNCTION_REVERSE_TWO_SPARSITY << "(j, &pos, &nnz);\n"
                "\n"
                "    compressed = (" << _baseTypeName << "*) malloc(nnz * sizeof (" << _baseTypeName << "));\n"
                "    in[0] = x;\n"
                "    in[1] = &tx[j * 2 + 1];\n"
                "    in[2] = py; \n"
                "    out[0] = compressed;\n"
                "    " << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO << "(j, in, out);\n"
                "\n"
                "    for (e = 0; e < nnz; e++) {\n"
                "        px[pos[e] * 2] = compressed[e];\n"
                "    }\n"
                "\n"
                "    free(compressed);\n"
                "    return 0;\n"
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
        _cache.str("");
        _cache << _name << "_" << function;
        std::string model_function = _cache.str();
        _cache.str("");
        generateFunctionDeclarationSource(_cache, model_function, suffix, elements);
        _cache << "\n";
        _cache << "int " << model_function << "("
                "unsigned long int pos," << _baseTypeName << " const *const * in, " << _baseTypeName << " *const * out) {\n"
                "   switch(pos) {\n";
        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            // the size of each sparsity row
            _cache << "      case " << it->first << ":\n"
                    "         " << model_function << "_" << suffix << it->first << "(in, out);\n"
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
    void CLangCompileModelHelper<Base>::generateFunctionDeclarationSource(std::ostringstream& cache,
                                                                          const std::string& model_function,
                                                                          const std::string& suffix,
                                                                          const std::map<size_t, std::vector<size_t> >& elements) {
        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            cache << "void " << model_function << "_" << suffix << it->first
                    << "(" << _baseTypeName << " const *const * in, " << _baseTypeName << " *const * out);\n";
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparsity1DSource(const std::string& function,
                                                                 const std::vector<size_t>& sparsity) {
        _cache << "void " << function << "("
                "unsigned long int const** sparsity,"
                " unsigned long int* nnz) {\n";

        // the size of each sparsity row
        _cache << "   static unsigned long int const nonzeros[" << sparsity.size() << "] = {";
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
                "unsigned long int const** row,"
                " unsigned long int const** col,"
                " unsigned long int* nnz) {\n";

        // the size of each sparsity row
        _cache << "static unsigned long int const rows[" << rows.size() << "] = {";
        if (!rows.empty()) {
            _cache << rows[0];
            for (size_t i = 1; i < rows.size(); i++) {
                _cache << "," << rows[i];
            }
        }
        _cache << "};\n";

        _cache << "static unsigned long int const cols[" << cols.size() << "] = {";
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
                "unsigned long int i,"
                "unsigned long int const** row,"
                " unsigned long int const** col,"
                " unsigned long int* nnz) {\n";

        for (size_t i = 0; i < sparsities.size(); i++) {
            const std::vector<size_t>& rows = sparsities[i].rows;
            const std::vector<size_t>& cols = sparsities[i].cols;
            assert(rows.size() == cols.size());

            _cache << "   static unsigned long int const rows" << i << "[" << rows.size() << "] = {";
            if (!rows.empty()) {
                _cache << rows[0];
                for (size_t i = 1; i < rows.size(); i++) {
                    _cache << "," << rows[i];
                }
            }
            _cache << "};\n";

            _cache << "   static unsigned long int const cols" << i << "[" << cols.size() << "] = {";
            if (!cols.empty()) {
                _cache << cols[0];
                for (size_t i = 1; i < cols.size(); i++) {
                    _cache << "," << cols[i];
                }
            }
            _cache << "};\n";
        }

        _cache << "   switch(i) {\n";
        for (size_t i = 0; i < sparsities.size(); i++) {
            // the size of each sparsity
            _cache << "   case " << i << ":\n"
                    "      *row = rows" << i << ";\n"
                    "      *col = cols" << i << ";\n"
                    "      *nnz = " << sparsities[i].rows.size() << ";\n"
                    "      break;\n";
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
                "unsigned long int pos,"
                " unsigned long int const** elements,"
                " unsigned long int* nnz) {\n";

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            // the size of each sparsity row
            const std::vector<size_t>& els = it->second;
            _cache << "   static unsigned long int const elements" << it->first << "[" << els.size() << "] = {";
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
            std::cout << "generating operation graph for '" << jobName << "' ... ";
            std::cout.flush();
            _beginTime = system::currentTime();
        }
    }

    template<class Base>
    void inline CLangCompileModelHelper<Base>::finishedGraphCreation() {
        if (_verbose) {
            double endTime = system::currentTime();
            std::cout << "done [" << std::fixed << std::setprecision(3)
                    << (endTime - _beginTime) << "]" << std::endl;
        }
    }

    template<class Base>
    inline std::map<size_t, std::vector<std::set<size_t> > > CLangCompileModelHelper<Base>::determineOrderByCol(const std::map<size_t, std::vector<size_t> >& elements,
                                                                                                                const LocalSparsityInfo& sparsity) {
        std::map<size_t, std::vector<std::set<size_t> > > userLocation;

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t col = it->first;
            const std::vector<size_t>& els = it->second;
            std::vector<std::set<size_t> >& userLocationCol = userLocation[col];
            userLocationCol.resize(els.size());

            for (size_t er = 0; er < els.size(); er++) {
                size_t row = els[er];
                for (size_t e = 0; e < sparsity.rows.size(); e++) {
                    if (sparsity.rows[e] == row && sparsity.cols[e] == col) {
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
        std::map<size_t, std::vector<std::set<size_t> > > userLocation;

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t row = it->first;
            const std::vector<size_t>& els = it->second;
            std::vector<std::set<size_t> >& userLocationRow = userLocation[row];
            userLocationRow.resize(els.size());

            for (size_t ec = 0; ec < els.size(); ec++) {
                size_t col = els[ec];
                for (size_t e = 0; e < sparsity.rows.size(); e++) {
                    if (sparsity.cols[e] == col && sparsity.rows[e] == row) {
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
