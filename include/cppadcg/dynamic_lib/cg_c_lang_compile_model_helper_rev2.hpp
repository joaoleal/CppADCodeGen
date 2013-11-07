#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_REV2_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_REV2_INCLUDED
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
    void CLangCompileModelHelper<Base>::generateSparseReverseTwoSources(std::map<std::string, std::string>& sources) throw (CGException) {

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

        if (!evalRows.empty()) {

            startingJob("'model (reverse two)'", JobTimer::SOURCE_GENERATION);

            if (isAtomicsUsed()) {
                generateSparseReverseTwoSourcesWithAtomics(sources, elements);
            } else {
                generateSparseReverseTwoSourcesNoAtomics(sources, elements, evalRows, evalCols);
            }

            finishedJob();

            _cache.str("");
        }

        generateGlobalDirectionalFunctionSource(FUNCTION_SPARSE_REVERSE_TWO,
                                                "indep",
                                                FUNCTION_REVERSE_TWO_SPARSITY,
                                                elements,
                                                sources);
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseReverseTwoSourcesWithAtomics(std::map<std::string, std::string>& sources,
                                                                                   const std::map<size_t, std::vector<size_t> >& elements) throw (CGException) {

        const size_t m = _fun.Range();
        const size_t n = _fun.Domain();
        //const size_t k = 1;
        const size_t p = 2;

        vector<CGBase> tx1v(n);

        std::map<size_t, std::vector<size_t> >::const_iterator it;
        for (it = elements.begin(); it != elements.end(); ++it) {
            size_t j = it->first;
            const std::vector<size_t>& cols = it->second;

            _cache.str("");
            _cache << "model (reverse two, indep " << j << ")";
            const std::string subJobName = _cache.str();

            startingJob("'" + subJobName + "'", JobTimer::GRAPH);

            CodeHandler<Base> handler;
            handler.setJobTimer(_jobTimer);

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

            finishedJob();

            CLanguage<Base> langC(_baseTypeName);
            langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
            langC.setParameterPrecision(_parameterPrecision);
            _cache.str("");
            _cache << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO << "_indep" << j;
            langC.setGenerateFunction(_cache.str());

            std::ostringstream code;
            std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("px"));
            CLangDefaultReverse2VarNameGenerator<Base> nameGenRev2(nameGen.get(), n, 1);

            handler.generateCode(code, langC, pxCustom, nameGenRev2, _atomicFunctions, subJobName);
        }
    }

    template<class Base>
    void CLangCompileModelHelper<Base>::generateSparseReverseTwoSourcesNoAtomics(std::map<std::string, std::string>& sources,
                                                                                 const std::map<size_t, std::vector<size_t> >& elements,
                                                                                 const std::vector<size_t>& evalRows,
                                                                                 const std::vector<size_t>& evalCols) throw (CGException) {
        const size_t m = _fun.Range();
        const size_t n = _fun.Domain();

        // save compressed positions
        std::map<size_t, std::map<size_t, size_t> > positions;
        std::map<size_t, std::vector<size_t> >::const_iterator itJ1;
        for (itJ1 = elements.begin(); itJ1 != elements.end(); ++itJ1) {
            size_t j1 = itJ1->first;
            const std::vector<size_t>& row = itJ1->second;
            std::map<size_t, size_t>& pos = positions[j1];

            for (size_t e = 0; e < row.size(); e++) {
                size_t j2 = row[e];
                pos[j2] = e;
            }
        }

        // we can use a new handler to reduce memory usage
        CodeHandler<Base> handler;
        handler.setJobTimer(_jobTimer);

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

        vector<CGBase> hessFlat(evalRows.size());

        CppAD::sparse_hessian_work work; // temporary structure for CPPAD
        _fun.SparseHessian(tx0, py, _hessSparsity.sparsity, evalRows, evalCols, hessFlat, work);

        std::map<size_t, vector<CGBase> > hess;
        for (itJ1 = elements.begin(); itJ1 != elements.end(); ++itJ1) {
            size_t j1 = itJ1->first;
            hess[j1].resize(itJ1->second.size());
        }

        // organize hessian elements
        for (size_t el = 0; el < evalRows.size(); el++) {
            size_t j1 = evalRows[el];
            size_t j2 = evalCols[el];
            size_t e = positions[j1][j2];

            hess[j1][e] = hessFlat[el];
        }

        /**
         * Generate one function for each independent variable
         */
        typename std::map<size_t, vector<CGBase> >::const_iterator it;
        for (it = hess.begin(); it != hess.end(); ++it) {
            size_t j = it->first;
            const vector<CGBase>& row = it->second;

            _cache.str("");
            _cache << "model (reverse two, indep " << j << ")";
            const std::string subJobName = _cache.str();

            vector<CGBase> pxCustom(row.size());
            for (size_t e = 0; e < row.size(); e++) {
                pxCustom[e] = row[e] * tx1;
            }

            CLanguage<Base> langC(_baseTypeName);
            langC.setMaxAssigmentsPerFunction(_maxAssignPerFunc, &sources);
            langC.setParameterPrecision(_parameterPrecision);
            _cache.str("");
            _cache << _name << "_" << FUNCTION_SPARSE_REVERSE_TWO << "_indep" << j;
            langC.setGenerateFunction(_cache.str());

            std::ostringstream code;
            std::auto_ptr<VariableNameGenerator<Base> > nameGen(createVariableNameGenerator("px"));
            CLangDefaultReverse2VarNameGenerator<Base> nameGenRev2(nameGen.get(), n, 1);

            handler.generateCode(code, langC, pxCustom, nameGenRev2, _atomicFunctions, subJobName);
        }
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

}

#endif
