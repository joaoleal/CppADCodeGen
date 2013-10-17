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
            _zeroEvaluated = true;
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

        startingJob("Loop detection");

        CodeHandler<Base> handler;
        handler.setJobTimer(this);

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

        finishedJob();
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
    bool CLangCompileModelHelper<Base>::isAtomicsUsed() {
        if (_zeroEvaluated) {
            return _atomicFunctions.size() > 0;
        } else {
            return !getAtomicsIndeps().empty();
        }
    }

    template<class Base>
    const std::map<size_t, std::set<size_t> >& CLangCompileModelHelper<Base>::getAtomicsIndeps() {
        if (_atomicsIndeps == NULL) {
            AtomicDependencyLocator<Base> adl(_fun);
            _atomicsIndeps = new std::map<size_t, std::set<size_t> >(adl.findAtomicsUsage());
        }
        return *_atomicsIndeps;
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
        _cache << "   ";
        CLanguage<Base>::printStaticIndexArray(_cache, "nonzeros", sparsity);

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
        _cache << "   ";
        CLanguage<Base>::printStaticIndexArray(_cache, "rows", rows);

        _cache << "   ";
        CLanguage<Base>::printStaticIndexArray(_cache, "cols", cols);

        _cache << "   *row = rows;\n"
                "   *col = cols;\n"
                "   *nnz = " << rows.size() << ";\n"
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

        std::ostringstream os;

        for (size_t i = 0; i < sparsities.size(); i++) {
            const std::vector<size_t>& rows = sparsities[i].rows;
            const std::vector<size_t>& cols = sparsities[i].cols;
            assert(rows.size() == cols.size());
            if (!rows.empty()) {
                os.str("");
                os << "rows" << i;
                _cache << "   ";
                CLanguage<Base>::printStaticIndexArray(_cache, os.str(), rows);

                os.str("");
                os << "cols" << i;
                _cache << "   ";
                CLanguage<Base>::printStaticIndexArray(_cache, os.str(), cols);
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
            _cache << "   ";
            std::ostringstream os;
            os << "elements" << it->first;
            CLanguage<Base>::printStaticIndexArray(_cache, os.str(), els);
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
