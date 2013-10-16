#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_FOR1_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_FOR1_INCLUDED
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
    void CLangCompileModelHelper<Base>::generateSparseForwardOneSources(std::map<std::string, std::string>& sources) {

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
        startingJob("source for 'model (reverse one)'");

        if (isAtomicsUsed()) {
            generateSparseForwardOneSourcesWithAtomics(sources, elements);
        } else {
            generateSparseForwardOneSourcesNoAtomics(sources, elements);
        }

        finishedJob();

        _cache.str("");

        generateGlobalDirectionalFunctionSource(FUNCTION_SPARSE_FORWARD_ONE,
                                                "indep",
                                                FUNCTION_FORWARD_ONE_SPARSITY,
                                                elements,
                                                sources);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseForwardOneSourcesWithAtomics(std::map<std::string, std::string>& sources,
                                                                                   const std::map<size_t, std::vector<size_t> >& elements) {
        /**
         * Generate one function for each dependent variable
         */
        size_t m = _fun.Range();
        size_t n = _fun.Domain();

        vector<CGBase> dxv(n);

        const std::string jobName = "model (forward one)";
        startingJob("source for '" + jobName + "'");

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t j = it->first;
            const std::vector<size_t>& rows = it->second;

            _cache.str("");
            _cache << "model (forward one, indep " << j << ")";
            const std::string subJobName = _cache.str();

            startingJob("operation graph for '" + subJobName + "'");

            CodeHandler<Base> handler;
            handler.setJobTimer(this);

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

            finishedJob();

            CLanguage<Base> langC(_baseTypeName);
            langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
            _cache.str("");
            _cache << _name << "_" << FUNCTION_SPARSE_FORWARD_ONE << "_indep" << j;
            langC.setGenerateFunction(_cache.str());

            std::ostringstream code;
            std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("dy", "x", "var", "array"));
            CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), "dx", n);

            handler.generateCode(code, langC, dyCustom, nameGenHess, _atomicFunctions, subJobName);
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseForwardOneSourcesNoAtomics(std::map<std::string, std::string>& sources,
                                                                                 const std::map<size_t, std::vector<size_t> >& elements) {
        /**
         * jacobian
         */
        size_t n = _fun.Domain();

        CodeHandler<Base> handler;
        handler.setJobTimer(this);

        vector<CGBase> x(n);
        handler.makeVariables(x);
        if (_x.size() > 0) {
            for (size_t i = 0; i < n; i++) {
                x[i].setValue(_x[i]);
            }
        }

        CGBase dx;
        handler.makeVariable(dx);
        if (_x.size() > 0) {
            dx.setValue(Base(1.0));
        }

        vector<CGBase> jacFlat(_jacSparsity.rows.size());

        CppAD::sparse_jacobian_work work; // temporary structure for CPPAD
        _fun.SparseJacobianForward(x, _jacSparsity.sparsity, _jacSparsity.rows, _jacSparsity.cols, jacFlat, work);

        /**
         * organize results
         */
        std::map<size_t, vector<CGBase> > jac; // by column
        std::map<size_t, std::map<size_t, size_t> > positions; // by column

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t j = it->first;
            const std::vector<size_t>& column = it->second;

            jac[j].resize(column.size());
            std::map<size_t, size_t>& pos = positions[j];

            for (size_t e = 0; e < column.size(); e++) {
                size_t i = column[e];
                pos[i] = e;
            }
        }

        for (size_t el = 0; el < _jacSparsity.rows.size(); el++) {
            size_t i = _jacSparsity.rows[el];
            size_t j = _jacSparsity.cols[el];
            size_t e = positions[j].at(i);

            vector<CGBase>& column = jac[j];
            column[e] = jacFlat[el] * dx;
        }

        /**
         * Create source for each independent/column
         */
        typename std::map<size_t, vector<CGBase> >::iterator itJ;
        for (itJ = jac.begin(); itJ != jac.end(); ++itJ) {
            size_t j = itJ->first;
            vector<CGBase>& dyCustom = itJ->second;

            _cache.str("");
            _cache << "model (forward one, indep " << j << ")";
            const std::string subJobName = _cache.str();

            CLanguage<Base> langC(_baseTypeName);
            langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
            _cache.str("");
            _cache << _name << "_" << FUNCTION_SPARSE_FORWARD_ONE << "_indep" << j;
            langC.setGenerateFunction(_cache.str());

            std::ostringstream code;
            std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("dy", "x", "var", "array"));
            CLangDefaultHessianVarNameGenerator<Base> nameGenHess(nameGen.get(), "dx", n);

            handler.generateCode(code, langC, dyCustom, nameGenHess, _atomicFunctions, subJobName);
        }
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
}

#endif
