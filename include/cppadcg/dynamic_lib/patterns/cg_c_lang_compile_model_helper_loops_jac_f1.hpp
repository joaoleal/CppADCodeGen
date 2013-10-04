#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_JAC_F1_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_JAC_F1_INCLUDED
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
    void CLangCompileModelHelper<Base>::generateSparseJacobianWithLoopsSourceFromFor1(std::map<std::string, std::string>& sources,
                                                                                      const std::map<size_t, std::vector<std::set<size_t> > >& userJacElLocation,
                                                                                      const std::map<size_t, bool>& jcolOrdered,
                                                                                      size_t maxCompressedSize) throw (CGException) {
        using namespace std;
        using namespace CppAD::loops;

        /**
         * Generate the source code
         */
        CLanguage<Base> langC(_baseTypeName);
        string argsDcl = langC.generateDefaultFunctionArgumentsDcl();

        string model_function = _name + "_" + FUNCTION_SPARSE_JACOBIAN;
        string functionFor1 = _name + "_" + FUNCTION_SPARSE_FORWARD_ONE;
        string suffix = "indep";
        string nlFor1Suffix = "noloop_" + suffix;
        IndexDclrOperationNode<Base> indexIt("it");

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n\n";
        
        generateFunctionDeclarationSource(_cache, functionFor1, nlFor1Suffix, _nonLoopFor1Elements, argsDcl);
        generateFunctionDeclarationSourceLoopForRev(_cache, langC, _name, "jcol", _loopFor1Groups, generateFunctionNameLoopFor1);
        
        _cache << "\n";
        printForRevUsageFunction<Base>(_cache, _baseTypeName, _name,
                model_function, 2,
                functionFor1, suffix,
                "jcol", indexIt, "jac",
                _loopFor1Groups,
                _nonLoopFor1Elements,
                userJacElLocation, jcolOrdered,
                CLangCompileModelHelper<Base>::generateFunctionNameLoopFor1,
                _jacSparsity.rows.size(), maxCompressedSize);

        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

}

#endif