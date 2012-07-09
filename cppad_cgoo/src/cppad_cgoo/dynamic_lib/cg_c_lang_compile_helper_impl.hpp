#ifndef CPPAD_CG_C_LANG_COMPILE_HELPER_IMPL_INCLUDED
#define	CPPAD_CG_C_LANG_COMPILE_HELPER_IMPL_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <typeinfo>

namespace CppAD {

    template<class Base>
    const unsigned long int CLangCompileHelper<Base>::API_VERSION = 0;

    template<class Base>
    const std::string CLangCompileHelper<Base>::FUNCTION_FORWAD_ZERO = "cppad_cg_forward_zero";

    template<class Base>
    const std::string CLangCompileHelper<Base>::FUNCTION_JACOBIAN = "cppad_cg_jacobian";

    template<class Base>
    const std::string CLangCompileHelper<Base>::FUNCTION_HESSIAN = "cppad_cg_hessian";

    template<class Base>
    const std::string CLangCompileHelper<Base>::FUNCTION_SPARSE_JACOBIAN = "cppad_cg_sparse_jacobian";

    template<class Base>
    const std::string CLangCompileHelper<Base>::FUNCTION_SPARSE_HESSIAN = "cppad_cg_sparse_hessian";

    template<class Base>
    const std::string CLangCompileHelper<Base>::FUNCTION_JACOBIAN_SPARSITY = "cppad_cg_jacobian_sparsity";

    template<class Base>
    const std::string CLangCompileHelper<Base>::FUNCTION_HESSIAN_SPARSITY = "cppad_cg_hessian_sparsity";

    template<class Base>
    const std::string CLangCompileHelper<Base>::FUNCTION_INFO = "cppad_cg_info";

    template<class Base>
    const std::string CLangCompileHelper<Base>::FUNCTION_VERSION = "cppad_cg_version";

    template<class Base>
    DynamicLib<Base>* CLangCompileHelper<Base>::createDynamicLibrary(CLangCompiler<Base>& compiler) {
        std::map<std::string, std::string> sources;
        if (_zero) {
            generateZeroSource();
            sources[FUNCTION_FORWAD_ZERO + ".c"] = _cache.str();
            _cache.str("");
        }

        if (_jacobian) {
            generateJacobianSource();
            sources[FUNCTION_JACOBIAN + ".c"] = _cache.str();
            _cache.str("");
        }

        if (_hessian) {
            generateHessianSource();
            sources[FUNCTION_HESSIAN + ".c"] = _cache.str();
            _cache.str("");
        }

        if (_sparseJacobian) {
            generateSparseJacobianSource(sources);
        }

        if (_sparseHessian) {
            generateSparseHessianSource(sources);
        }

        generateVerionSource();
        sources[FUNCTION_VERSION + ".c"] = _cache.str();
        _cache.str("");

        generateInfoSource();
        sources[FUNCTION_INFO + ".c"] = _cache.str();
        _cache.str("");

        compiler.compileDynamic(_libraryName, sources, true);

        return loadDynamicLibrary();
    }

    template<class Base>
    void CLangCompileHelper<Base>::generateVerionSource() {
        _cache << "unsigned long int " << FUNCTION_VERSION << "() {\n";
        _cache << "return " << API_VERSION << "u;\n";
        _cache << "}\n\n";
    }

    template<class Base>
    void CLangCompileHelper<Base>::generateInfoSource() {
        const char* localBaseName = typeid (Base).name();

        _cache << "void " << FUNCTION_INFO << "(const char** baseName, unsigned long int* m, unsigned long int* n) {\n";
        _cache << "*baseName = \"" << baseTypeName() << "  " << localBaseName << "\";\n";
        _cache << "*m = " << _fun->Range() << ";\n";
        _cache << "*n = " << _fun->Domain() << ";\n";
        _cache << "}\n\n";
    }

    template<class Base>
    void CLangCompileHelper<Base>::generateZeroSource() {
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCG;

        CodeHandler<Base> handler;

        std::vector<CGD> indVars(_fun->Domain());
        handler.makeVariables(indVars);

        std::vector<CGD> dep = _fun->Forward(0, indVars);

        CLanguage<Base> langC;
        CLangDefaultVariableNameGenerator<Base> nameGen;

        std::ostringstream code;
        handler.generateCode(code, langC, dep, nameGen);

        std::string basename = baseTypeName();

        _cache << "#include <math.h>\n\n";
        _cache << "void " << FUNCTION_FORWAD_ZERO << "(const " << basename << "* ind, " << basename << "* dep) {\n";
        // declare temporary variables
        _cache << langC.generateTemporaryVariableDeclaration(basename);
        // the code
        _cache << code.str();
        _cache << "}\n\n";
    }

    template<class Base>
    void CLangCompileHelper<Base>::generateJacobianSource() {
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCG;

        CodeHandler<Base> handler;

        std::vector<CGD> indVars(_fun->Domain());
        handler.makeVariables(indVars);

        std::vector<CGD> jac = _fun->Jacobian(indVars);

        CLanguage<Base> langC;
        CLangDefaultVariableNameGenerator<Base> nameGen("jac", "ind", "var");

        std::ostringstream code;
        handler.generateCode(code, langC, jac, nameGen);

        std::string basename = baseTypeName();

        _cache << "#include <math.h>\n\n";
        _cache << "void " << FUNCTION_JACOBIAN << "(const " << basename << "* ind, " << basename << "* jac) {\n";
        // declare temporary variables
        _cache << langC.generateTemporaryVariableDeclaration(basename);
        // the code
        _cache << code.str();
        _cache << "}\n\n";
    }

    template<class Base>
    void CLangCompileHelper<Base>::generateHessianSource() {
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCG;

        CodeHandler<Base> handler;

        size_t m = _fun->Domain();

        // independent variables
        std::vector<CGD> indVars(m);
        handler.makeVariables(indVars);
        // multipliers
        std::vector<CGD> w(_fun->Range());
        handler.makeVariables(w);

        std::vector<CGD> hess = _fun->Hessian(indVars, w);

        // make use of the symmetry of the Hessian in order to reduce operations
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < i; j++) {
                hess[i * m + j] = hess[j * m + i];
            }
        }

        CLanguage<Base> langC;
        CLangDefaultHessianVarNameGenerator<Base> nameGen(_fun->Domain());

        std::ostringstream code;
        handler.generateCode(code, langC, hess, nameGen);

        std::string basename = baseTypeName();

        _cache << "#include <math.h>\n\n";
        _cache << "void " << FUNCTION_HESSIAN << "(const " << basename << "* ind, const " << basename << "* mult, " << basename << "* hess) {\n";
        // declare temporary variables
        _cache << langC.generateTemporaryVariableDeclaration(basename);
        // the code
        _cache << code.str();
        _cache << "}\n\n";
    }

    template<class Base>
    void CLangCompileHelper<Base>::generateSparseJacobianSource(std::map<std::string, std::string>& sources) {
        size_t m = _fun->Range();
        size_t n = _fun->Domain();

        std::vector<size_t> rows, cols;

        /**
         * Determine the sparsity pattern
         */
        std::vector<bool> sparsity;
        if (n <= m) {
            // use forward mode 
            std::vector<bool> r(n * n);
            for (size_t j = 0; j < n; j++) {
                for (size_t k = 0; k < n; k++)
                    r[j * n + k] = false;
                r[j * n + j] = true;
            }
            sparsity = _fun->ForSparseJac(n, r);
        } else {
            // use reverse mode 
            std::vector<bool> s(m * m);
            for (size_t i = 0; i < m; i++) {
                for (size_t k = 0; k < m; k++)
                    s[i * m + k] = false;
                s[i * m + i] = true;
            }
            sparsity = _fun->RevSparseJac(m, s);
        }

        if (_custom_jac_row.empty()) {
            generateSparsityIndexes(sparsity, m, n, rows, cols);

        } else {
            rows = _custom_jac_row;
            cols = _custom_jac_col;
        }

        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCG;

        CodeHandler<Base> handler;

        std::vector<CGD> indVars(n);
        handler.makeVariables(indVars);

        std::vector<CGD> jac(rows.size());
        CppAD::sparse_jacobian_work work;
        if (n <= m) {
            _fun->SparseJacobianForward(indVars, sparsity, rows, cols, jac, work);
        } else {
            _fun->SparseJacobianReverse(indVars, sparsity, rows, cols, jac, work);
        }

        CLanguage<Base> langC;
        CLangDefaultVariableNameGenerator<Base> nameGen("jac", "ind", "var");

        std::ostringstream code;
        handler.generateCode(code, langC, jac, nameGen);

        std::string basename = baseTypeName();

        _cache << "#include <math.h>\n\n";
        _cache << "void " << FUNCTION_SPARSE_JACOBIAN << "(const " << basename << "* ind, " << basename << "* jac) {\n";
        // declare temporary variables
        _cache << langC.generateTemporaryVariableDeclaration(basename);
        // the code
        _cache << code.str();
        _cache << "}\n\n";

        sources[FUNCTION_SPARSE_JACOBIAN + ".c"] = _cache.str();
        _cache.str("");

        generateSparsitySource(FUNCTION_JACOBIAN_SPARSITY, rows, cols);
        sources[FUNCTION_JACOBIAN_SPARSITY + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileHelper<Base>::generateSparseHessianSource(std::map<std::string, std::string>& sources) {
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

        /**
         * 
         */
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCG;

        CodeHandler<Base> handler;

        // independent variables
        std::vector<CGD> indVars(n);
        handler.makeVariables(indVars);
        // multipliers
        std::vector<CGD> w(m);
        handler.makeVariables(w);

        CppAD::sparse_hessian_work work;
        std::vector<CGD> hess(rows.size());
        _fun->SparseHessian(indVars, w, sparsity, rows, cols, hess, work);

        CLanguage<Base> langC;
        CLangDefaultHessianVarNameGenerator<Base> nameGen(n);

        std::ostringstream code;
        handler.generateCode(code, langC, hess, nameGen);

        std::string basename = baseTypeName();

        _cache << "#include <math.h>\n\n";
        _cache << "void " << FUNCTION_SPARSE_HESSIAN << "(const " << basename << "* ind, const " << basename << "* mult, " << basename << "* hess) {\n";
        // declare temporary variables
        _cache << langC.generateTemporaryVariableDeclaration(basename);
        // the code
        _cache << code.str();
        _cache << "}\n\n";

        sources[FUNCTION_SPARSE_HESSIAN + ".c"] = _cache.str();
        _cache.str("");

        generateSparsitySource(FUNCTION_HESSIAN_SPARSITY, rows, cols);
        sources[FUNCTION_HESSIAN_SPARSITY + ".c"] = _cache.str();
        _cache.str("");
    }

    template<class Base>
    void CLangCompileHelper<Base>::generateSparsitySource(const std::string& function,
                                                          const std::vector<bool>& sparsity,
                                                          size_t m, size_t n) {
        std::vector<size_t> rows, cols;

        generateSparsityIndexes(sparsity, m, n, rows, cols);

        generateSparsitySource(function, rows, cols);
    }

    template<class Base>
    void CLangCompileHelper<Base>::generateSparsitySource(const std::string& function,
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
    void CLangCompileHelper<Base>::generateSparsityIndexes(const std::vector<bool>& sparsity,
                                                           size_t m, size_t n,
                                                           std::vector<size_t>& rows,
                                                           std::vector<size_t>& cols) {
        size_t K = 0;
        for (size_t i = 0; i < sparsity.size(); i++) {
            if (sparsity[i]) {
                K++;
            }
        }

        rows.resize(K);
        cols.resize(K);

        size_t k = 0;
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                if (sparsity[i * n + j]) {
                    rows[k] = i;
                    cols[k] = j;
                    k++;
                }
            }
        }

        assert(k == K);
    }

    /**
     * 
     * Specializations
     */
    template<>
    inline const std::string CLangCompileHelper<double>::baseTypeName() {
        return "double";
    }

    template<>
    inline const std::string CLangCompileHelper<float>::baseTypeName() {
        return "float";
    }
}

#endif
