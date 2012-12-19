#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_IMPL_INCLUDED
#define	CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_IMPL_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

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
    const std::string CLangCompileModelHelper<Base>::FUNCTION_SPARSE_JACOBIAN = "sparse_jacobian";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_SPARSE_HESSIAN = "sparse_hessian";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_JACOBIAN_SPARSITY = "jacobian_sparsity";

    template<class Base>
    const std::string CLangCompileModelHelper<Base>::FUNCTION_HESSIAN_SPARSITY = "hessian_sparsity";

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
    void CLangCompileModelHelper<Base>::compileSources(CLangCompiler<Base>& compiler) {
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

        generateInfoSource(sources);

        compiler.compileSources(sources, true);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateInfoSource(std::map<std::string, std::string>& sources) {
        const char* localBaseName = typeid (Base).name();

        std::string funcName = _name + "_" + FUNCTION_INFO;

        std::auto_ptr<VariableNameGenerator< Base > > nameGen(createVariableNameGenerator("dep", "ind", "var"));

        _cache.str("");
        _cache << "void " << funcName << "(const char** baseName, unsigned long int* m, unsigned long int* n, unsigned int* indCount, unsigned int* depCount) {\n";
        _cache << "   *baseName = \"" << _baseTypeName << "  " << localBaseName << "\";\n";
        _cache << "   *m = " << _fun->Range() << ";\n";
        _cache << "   *n = " << _fun->Domain() << ";\n";
        _cache << "   *depCount = " << nameGen->getDependent().size() << ";\n";
        _cache << "   *indCount = " << nameGen->getIndependent().size() << ";\n";
        _cache << "}\n\n";

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
        const std::string jobName = "sparse Jacobian";

        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        std::vector<size_t> rows, cols;

        /**
         * Determine the sparsity pattern
         */
        std::vector<bool> sparsity = jacobianSparsity<std::vector<bool>, CGBase>(*_fun);

        if (_custom_jac_row.empty()) {
            generateSparsityIndexes(sparsity, m, n, rows, cols);

        } else {
            rows = _custom_jac_row;
            cols = _custom_jac_col;
        }

        startingGraphCreation(jobName);

        CodeHandler<Base> handler;
        handler.setVerbose(_verbose);

        std::vector<CGBase> indVars(n);
        handler.makeVariables(indVars);

        std::vector<CGBase> jac(rows.size());
        CppAD::sparse_jacobian_work work;
        if (n <= m) {
            _fun->SparseJacobianForward(indVars, sparsity, rows, cols, jac, work);
        } else {
            _fun->SparseJacobianReverse(indVars, sparsity, rows, cols, jac, work);
        }

        finishedGraphCreation();

        CLanguage<Base> langC(_baseTypeName);
        langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
        langC.setGenerateFunction(_name + "_" + FUNCTION_SPARSE_JACOBIAN);

        std::ostringstream code;
        std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("jac", "ind", "var"));

        handler.generateCode(code, langC, jac, *nameGen, jobName);

        generateSparsitySource(_name + "_" + FUNCTION_JACOBIAN_SPARSITY, rows, cols);
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
        std::vector<bool> s(m, true);
        std::vector<bool> r(n * n);
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < n; k++)
                r[j * n + k] = false;
            r[j * n + j] = true;
        }
        _fun->ForSparseJac(n, r);
        std::vector<bool> sparsity = _fun->RevSparseHes(n, s);

        std::vector<size_t> rows, cols;
        if (_custom_hess_row.empty()) {
            generateSparsityIndexes(sparsity, n, n, rows, cols);

        } else {
            rows = _custom_hess_row;
            cols = _custom_hess_col;
        }

        // make use of the symmetry of the Hessian in order to reduce operations
        std::map<size_t, std::map<size_t, size_t> > locations;
        for (size_t i = 0; i < rows.size(); i++) {
            locations[rows[i]][cols[i]] = i;
        }

        std::vector<size_t> upperHessRows, upperHessCols, upperHessOrder;
        upperHessRows.reserve(rows.size() / 2);
        upperHessCols.reserve(upperHessRows.size());
        upperHessOrder.reserve(upperHessRows.size());

        std::map<size_t, size_t> duplicates; // the elements determined using symmetry
        std::map<size_t, std::map<size_t, size_t> >::const_iterator ii;
        std::map<size_t, size_t>::const_iterator jj;
        for (size_t i = 0; i < rows.size(); i++) {
            bool add = true;
            if (rows[i] > cols[i]) {
                ii = locations.find(cols[i]);
                if (ii != locations.end()) {
                    jj = ii->second.find(cols[i]);
                    if (jj != ii->second.end()) {
                        size_t k = jj->second;
                        duplicates[i] = k;
                        add = false; // symmetric value being determined
                    }
                }
            }

            if (add) {
                upperHessRows.push_back(rows[i]);
                upperHessCols.push_back(cols[i]);
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
        _fun->SparseHessian(indVars, w, sparsity, upperHessRows, upperHessCols, upperHess, work);

        std::vector<CGBase> hess(rows.size());
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

        generateSparsitySource(_name + "_" + FUNCTION_HESSIAN_SPARSITY, rows, cols);
        sources[_name + "_" + FUNCTION_HESSIAN_SPARSITY + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparsitySource(const std::string& function,
                                                               const std::vector<size_t>& rows,
                                                               const std::vector<size_t>& cols) {
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
